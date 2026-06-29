#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import sys
import json
import warnings

from rdkit import Chem

from meeko import MoleculePreparation
from meeko import rdkitutils
from meeko.pdbt_writer import PDBTWriterLegacy
from meeko.pdbt_atomtyper import assign_pdbt_types


def cmd_lineparser():
    conf_parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
    )
    conf_parser.add_argument(
        "-c", "--config_file",
        help="configure MoleculePreparation from JSON file. Overriden by command line args.",
    )
    confargs, remaining_argv = conf_parser.parse_known_args()

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-v", "--verbose",
        dest="verbose", action="store_true",
        help="print information about molecule setup",
    )

    io_group = parser.add_argument_group("Input/Output")
    io_group.add_argument(
        "-i", "--mol",
        dest="input_molecule_filename", required=True, action="store",
        help="molecule file (MOL2, SDF, ...)",
    )
    io_group.add_argument(
        "-o", "--out",
        dest="output_pdbt_filename", action="store",
        help="output pdbt filename. Single molecule input only.",
    )
    io_group.add_argument(
        "-", "--",
        dest="redirect_stdout", action="store_true",
        help="do not write file, redirect output to STDOUT.",
    )

    config_group = parser.add_argument_group("Molecule preparation")
    config_group.add_argument(
        "-c", "--config_file",
        help="configure MoleculePreparation from JSON file.",
    )
    config_group.add_argument(
        "--rigid_macrocycles",
        dest="rigid_macrocycles", action="store_true",
        help="keep macrocycles rigid in input conformation",
    )
    config_group.add_argument(
        "--keep_chorded_rings",
        dest="keep_chorded_rings", action="store_true",
        help="return all rings from exhaustive perception",
    )
    config_group.add_argument(
        "--keep_equivalent_rings",
        dest="keep_equivalent_rings", action="store_true",
        help="equivalent rings have the same size and neighbors",
    )
    config_group.add_argument(
        "--min_ring_size",
        dest="min_ring_size", type=int,
        help="min nr of atoms in ring for opening",
    )
    config_group.add_argument(
        "-r", "--rigidify_bonds_smarts",
        dest="rigidify_bonds_smarts", action="append",
        help="SMARTS patterns to rigidify bonds", metavar="SMARTS",
    )
    config_group.add_argument(
        "-b", "--rigidify_bonds_indices",
        dest="rigidify_bonds_indices", action="append",
        help="indices of two atoms (in the SMARTS) that define a bond (start at 1)",
        nargs="+", type=int, metavar="i j",
    )
    config_group.add_argument(
        "-a", "--flexible_amides",
        dest="flexible_amides", action="store_true",
        help="allow amide bonds to rotate and be non-planar, which is bad",
    )
    config_group.add_argument(
        "--double_bond_penalty",
        help="penalty > 100 prevents breaking double bonds",
        type=int,
    )
    config_group.add_argument(
        "--charge_model",
        choices=("gasteiger", "zero", "read"),
        help="default is 'gasteiger', 'zero' sets all zeros",
        default="gasteiger",
    )
    config_group.add_argument(
        "--charge_atom_prop",
        help="set atom partial charges from an RDKit atom property based on the input file.",
    )
    config_group.add_argument(
        "--bad_charge_ok",
        help="NaN and Inf charges allowed in PDBT",
        action="store_true",
    )
    config_group.add_argument(
        "--rename_atoms",
        dest="rename_atoms", action="store_true",
        help="rename atoms: new name is original name + (1-based) index",
    )

    config = MoleculePreparation.get_defaults_dict()
    if confargs.config_file is not None:
        with open(confargs.config_file) as f:
            c = json.load(f)
            config.update(c)

    parser.set_defaults(**config)
    args = parser.parse_args(remaining_argv)

    for key in config:
        if key in args.__dict__:
            config[key] = args.__dict__[key]

    return args, config


def main():
    args, config = cmd_lineparser()
    input_molecule_filename = args.input_molecule_filename

    input_fname, ext = os.path.splitext(input_molecule_filename)
    ext = ext[1:].lower()

    parsers = {
        "sdf": Chem.SDMolSupplier,
        "mol2": rdkitutils.Mol2MolSupplier,
        "mol": Chem.SDMolSupplier,
    }
    if ext not in parsers:
        print(
            "Error: Format [%s] not in supported formats [%s]"
            % (ext, "/".join(list(parsers.keys())))
        )
        sys.exit(1)
    mol_supplier = parsers[ext](input_molecule_filename, removeHs=False)

    if args.output_pdbt_filename is None:
        output_filename = input_fname + ".pdbt"
    else:
        output_filename = args.output_pdbt_filename

    if args.charge_atom_prop is not None:
        if config["charge_model"] != "read":
            print(
                'Error: --charge_atom_prop must be used with --charge_model "read"',
                file=sys.stderr,
            )
            sys.exit(1)
    elif config["charge_model"] == "read":
        if ext == "sdf":
            config["charge_atom_prop"] = "PartialCharge"
        elif ext == "mol2":
            config["charge_atom_prop"] = "_TriposPartialCharge"

    preparator = MoleculePreparation.from_config(config)

    nr_failures = 0
    is_after_first = False

    for mol in mol_supplier:
        if is_after_first:
            print("Processed only the first molecule of multiple molecule input.", file=sys.stderr)
            break

        if mol is None:
            continue

        name = mol.GetProp("_Name")
        is_after_first = True

        try:
            molsetups = preparator.prepare(mol, rename_atoms=args.rename_atoms)
        except Exception as error_msg:
            nr_failures += 1
            print(error_msg, file=sys.stderr)
            continue

        for molsetup in molsetups:
            assign_pdbt_types(molsetup, mol)

            pdbt_string, success, error_msg = PDBTWriterLegacy.write_string(
                molsetup,
                bad_charge_ok=args.bad_charge_ok,
            )

            if success:
                if args.redirect_stdout:
                    print(pdbt_string, end="")
                else:
                    out_fn = output_filename if len(molsetups) == 1 else "%s_%s.pdbt" % (input_fname, molsetup.name)
                    with open(out_fn, "w") as f:
                        f.write(pdbt_string)
                    print("Wrote %s" % out_fn)
                if args.verbose:
                    molsetup.show()
            else:
                nr_failures += 1
                print(error_msg, file=sys.stderr)

    if nr_failures > 0:
        sys.exit(1)


if __name__ == "__main__":
    sys.exit(main())
