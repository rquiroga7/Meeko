#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Prepare a receptor in PDBT (Vinardo2) format for AutoDock docking.

This is a thin, focused CLI for writing a single rigid-body receptor
PDBT. It is independent from ``mk_prepare_receptor.py``: no flexres,
no reactive residues, no GPF/box output. The idea is to mirror
``meeko-pdbt`` (the ligand CLI) for the receptor side.

Examples
--------

Prepare a PDB receptor and write a PDBT file::

    meeko-pdbt-rec -i receptor.pdb -o receptor.pdbt

Read a mmCIF file (requires ProDy) and write to stdout::

    meeko-pdbt-rec -i receptor.cif -o -

Handle alternate locations::

    meeko-pdbt-rec -i receptor.pdb -o receptor.pdbt --default_altloc A
"""

import argparse
import json
import os
import sys
import warnings
import pathlib

from meeko import MoleculePreparation
from meeko import PDBTWriterLegacy
from meeko import Polymer
from meeko import PolymerCreationError
from meeko import ResidueChemTemplates
from meeko.pdbt_atomtyper import assign_pdbt_types_from_pdbinfo

try:
    import prody
except ImportError as import_error:
    _prody_import_error = import_error
    _got_prody = False
else:
    SUPPORTED_PRODY_FORMATS = {"pdb": prody.parsePDB, "cif": prody.parseMMCIF}
    _got_prody = True


def cmd_lineparser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    io_group = parser.add_argument_group("Input/Output")
    io_group.add_argument(
        "-i", "--input",
        dest="input_filename", required=True,
        help="receptor file (PDB or mmCIF; mmCIF requires ProDy)",
    )
    io_group.add_argument(
        "-o", "--output",
        dest="output_filename", default=None,
        help="output PDBT filename. Default: <input>.pdbt. Use '-' for stdout.",
    )
    io_group.add_argument(
        "-", "--stdout",
        dest="redirect_stdout", action="store_true",
        help="do not write file, redirect output to STDOUT",
    )

    config_group = parser.add_argument_group("Receptor perception")
    config_group.add_argument(
        "--default_altloc",
        help="default alternate location (e.g. 'A')",
    )
    config_group.add_argument(
        "--wanted_altloc",
        help="require altloc for specific residues, e.g. :5=B,B:17=A",
    )
    config_group.add_argument(
        "-n", "--set_template",
        help="override residue template, e.g. A:5,7=CYX",
    )
    config_group.add_argument(
        "-d", "--delete_residues",
        help="delete residues by chain:number, e.g. A:350,B:15,16,17",
    )
    config_group.add_argument(
        "-b", "--blunt_ends",
        help="mark chain ends as blunt, e.g. A:123=2,A:1=0",
    )
    config_group.add_argument(
        "--config_file",
        help="JSON file with extra MoleculePreparation settings",
    )
    config_group.add_argument(
        "--add_templates",
        action="append", default=[],
        help="additional residue template (JSON file or 'resname:file.sdf'); repeatable",
    )
    config_group.add_argument(
        "-x", "--delete_bad_res",
        action="store_true",
        help="delete residues that don't match templates instead of raising",
    )
    config_group.add_argument(
        "--forgive_extra_bonds",
        action="store_true",
        help="allow processing structures with excess bonds (use with care)",
    )
    config_group.add_argument(
        "--bad_charge_ok",
        action="store_true",
        help="allow NaN/Inf charges in the output",
    )
    return parser.parse_args()


def _parse_blunt_ends(value):
    if value is None:
        return []
    out = []
    for chunk in value.split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        if "=" not in chunk:
            raise ValueError(f"bad --blunt_ends entry {chunk!r}; expected resid=value")
        resid, val = chunk.split("=", 1)
        out.append((resid.strip(), int(val.strip())))
    return out


def _parse_resid_list(value):
    """Parse 'A:5,B:7' or 'A:5,6,17' into {'A:5', 'A:6', 'A:17'}."""
    out = set()
    for chunk in value.split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        if ":" in chunk:
            chain, nums = chunk.split(":", 1)
            chain = chain.strip()
            for n in nums.split(","):
                n = n.strip()
                if n:
                    out.add(f"{chain}:{n}")
        else:
            out.add(f":{chunk}")
    return out


def _parse_set_template(value):
    """Parse 'A:5=CYX,B:17=HID' into {'A:5': 'CYX', 'B:17': 'HID'}."""
    out = {}
    for chunk in value.split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        if "=" not in chunk:
            raise ValueError(f"bad --set_template entry {chunk!r}; expected resid=resname")
        resid, name = chunk.split("=", 1)
        out[resid.strip()] = name.strip()
    return out


def _build_polymer(args):
    """Build a Polymer from the input file using the appropriate reader."""
    input_path = pathlib.Path(args.input_filename)
    if not input_path.exists():
        raise FileNotFoundError(input_path)
    suffix = input_path.suffix.lower()

    set_template = _parse_set_template(args.set_template) if args.set_template else {}
    delete_residues = _parse_resid_list(args.delete_residues) if args.delete_residues else []
    blunt_ends = _parse_blunt_ends(args.blunt_ends) if args.blunt_ends else []
    wanted_altloc = (
        {r.strip(): v.strip() for chunk in args.wanted_altloc.split(",") if "=" in chunk
         for r, v in [chunk.split("=", 1)]}
        if args.wanted_altloc else None
    )

    # Build mk_prep
    if args.config_file:
        with open(args.config_file) as f:
            config = json.load(f)
    else:
        config = {}
    mk_prep = MoleculePreparation.from_config(config) if config else MoleculePreparation()

    # Build templates (with optional --add_templates)
    templates = ResidueChemTemplates.create_from_defaults()
    for item in args.add_templates:
        if item.endswith(".json"):
            templates.add_json_file(item)
        elif ":" in item:
            from rdkit import Chem as _Chem
            from meeko.cli.mk_prepare_receptor import sdf_to_json
            resname, sdf_file = item.split(":", 1)
            templates.add_dict(sdf_to_json(sdf_file, resname))
        else:
            raise ValueError(f"--add_templates entry must be a JSON file or 'resname:file.sdf', got {item!r}")

    if suffix in (".cif", ".mmcif"):
        if not _got_prody:
            raise RuntimeError(
                f"reading {suffix} requires ProDy. Install with "
                "`pip install prody` or convert the input to PDB."
            )
        prody_obj = SUPPORTED_PRODY_FORMATS["cif"](str(input_path), altloc="all")
        return Polymer.from_prody(
            prody_obj, templates, mk_prep, set_template,
            delete_residues, False, args.delete_bad_res,
            blunt_ends=blunt_ends,
            wanted_altloc=wanted_altloc,
            default_altloc=args.default_altloc,
            forgive_extra_bonds=args.forgive_extra_bonds,
        )
    elif suffix == ".pdb":
        with open(input_path) as f:
            pdb_string = f.read()
        return Polymer.from_pdb_string(
            pdb_string, templates, mk_prep, set_template,
            delete_residues, False, args.delete_bad_res,
            blunt_ends=blunt_ends,
            wanted_altloc=wanted_altloc,
            default_altloc=args.default_altloc,
            forgive_extra_bonds=args.forgive_extra_bonds,
        )
    else:
        raise ValueError(
            f"unsupported input format {suffix!r}; use .pdb or .cif/.mmcif"
        )


def main():
    args = cmd_lineparser()

    try:
        polymer = _build_polymer(args)
    except PolymerCreationError as e:
        print(e, file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError as e:
        print(f"input file not found: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"error reading receptor: {e}", file=sys.stderr)
        sys.exit(1)

    n_monomers = len(polymer.get_valid_monomers())
    print(f"Read {n_monomers} residues from {args.input_filename}", file=sys.stderr)

    # Assign PDBT types using the hardcoded protein lookup
    assign_pdbt_types_from_pdbinfo(polymer)

    # Write the PDBT
    pdbt_string, success, error_msg = PDBTWriterLegacy.write_from_polymer(
        polymer, bad_charge_ok=args.bad_charge_ok,
    )
    if not success:
        print(error_msg, file=sys.stderr)
        sys.exit(1)

    # Decide where to write
    if args.redirect_stdout or args.output_filename == "-":
        sys.stdout.write(pdbt_string)
    else:
        out_name = args.output_filename
        if out_name is None:
            out_name = str(pathlib.Path(args.input_filename).with_suffix(".pdbt"))
        with open(out_name, "w") as f:
            f.write(pdbt_string)
        print(f"Wrote {out_name}", file=sys.stderr)


if __name__ == "__main__":
    sys.exit(main())
