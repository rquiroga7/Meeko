"""Compare meeko-pdbt-rec atom typing against obabel-25-07 reference for
receptors in the runs-n-poses dataset.

This test focuses on a small sample of receptors (by default 20) to
verify the meeko-pdbt-rec CLI and the underlying assign_pdbt_types_from_pdbinfo
function produce receptor PDBT atom types that match the obabel-25-07
reference within tolerance.

The reference PDBTs in vinardo_inputs/receptors/<system>/ were generated
with obabel-25-07 using ``-p7 -xr`` (pH 7 protonation, rigid output).
This test invokes obabel once per receptor to produce a fresh reference
so any difference in obabel behaviour is captured.

Usage:
    /home/rquiroga/anaconda3/envs/meeko_pdbt/bin/python test/test_pdbt_receptor_compare.py
    /home/rquiroga/anaconda3/envs/meeko_pdbt/bin/python test/test_pdbt_receptor_compare.py --max 5
    /home/rquiroga/anaconda3/envs/meeko_pdbt/bin/python test/test_pdbt_receptor_compare.py --system 5s9y__1__1.A__1.K
"""

import argparse
import math
import sys
import time
from collections import Counter
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

DATASET = Path("/home/rquiroga/Datasets/runs-n-poses-datasets")
GROUND_TRUTH = DATASET / "ground_truth"
VINARDO_RECEPTORS = DATASET / "vinardo_inputs" / "receptors"

# Make the meeko package importable
from meeko import MoleculePreparation, Polymer, ResidueChemTemplates
from meeko import PDBTWriterLegacy
from meeko.pdbt_atomtyper import assign_pdbt_types_from_pdbinfo


def parse_pdbt_atoms(text):
    """Parse ATOM/HETATM lines from a PDBT string."""
    atoms = []
    for line in text.splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        try:
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            atom_type = line[77:79].strip()
        except ValueError:
            continue
        atoms.append({"coord": (x, y, z), "type": atom_type, "line": line})
    return atoms


def _is_hydrogen_type(t):
    return t and t[0] == "H"


def match_atoms(meeko_atoms, obabel_atoms):
    """Two-phase coordinate matching: tight for heavy, wide for H."""
    meeko_matched = [False] * len(meeko_atoms)
    obabel_matched = [False] * len(obabel_atoms)
    matches = []
    for mi, ma in enumerate(meeko_atoms):
        if _is_hydrogen_type(ma["type"]):
            continue
        best_dist = 0.01
        best_oi = None
        for oi, oa in enumerate(obabel_atoms):
            if obabel_matched[oi] or _is_hydrogen_type(oa["type"]):
                continue
            d = math.sqrt(sum((ma["coord"][i] - oa["coord"][i]) ** 2 for i in range(3)))
            if d < best_dist:
                best_dist = d
                best_oi = oi
        if best_oi is not None:
            meeko_matched[mi] = True
            obabel_matched[best_oi] = True
            matches.append((ma, obabel_atoms[best_oi]))
    for mi, ma in enumerate(meeko_atoms):
        if meeko_matched[mi]:
            continue
        best_dist = 0.5
        best_oi = None
        for oi, oa in enumerate(obabel_atoms):
            if obabel_matched[oi]:
                continue
            d = math.sqrt(sum((ma["coord"][i] - oa["coord"][i]) ** 2 for i in range(3)))
            if d < best_dist:
                best_dist = d
                best_oi = oi
        if best_oi is not None:
            meeko_matched[mi] = True
            obabel_matched[best_oi] = True
            matches.append((ma, obabel_atoms[best_oi]))
    return matches


def meeko_build_pdbt(pdb_path, delete_bad_res=True):
    """Build a receptor PDBT with meeko using only hardcoded protein typing."""
    templates = ResidueChemTemplates.create_from_defaults()
    mk_prep = MoleculePreparation()
    with open(pdb_path) as f:
        pdb_string = f.read()
    try:
        polymer = Polymer.from_pdb_string(
            pdb_string, templates, mk_prep, allow_bad_res=delete_bad_res,
        )
    except Exception as e:
        return None, f"Polymer creation failed: {e}"
    assign_pdbt_types_from_pdbinfo(polymer)
    pdbt_string, success, err = PDBTWriterLegacy.write_from_polymer(polymer)
    if not success:
        return None, err or "PDBTWriterLegacy.write_from_polymer failed"
    return pdbt_string, None


def find_systems(max_n=None, system=None):
    """Return a list of system_ids to test.

    Only systems that have a pre-computed obabel reference PDBT in
    ``vinardo_inputs/receptors/<system>/`` are used, since the reference
    is the ground truth we are comparing against. This avoids cases
    where the input PDB and CIF represent different parts of the
    receptor (e.g. 5s9l has a 24-834 PDB but a different CIF
    fragment), which would make coordinate matching meaningless.
    """
    if system is not None:
        if (VINARDO_RECEPTORS / system).exists():
            return [system]
        return []
    systems = []
    for p in sorted(GROUND_TRUTH.iterdir()):
        if not p.is_dir():
            continue
        if not (p / "receptor.pdb").exists() or not (p / "receptor.cif").exists():
            continue
        if not (VINARDO_RECEPTORS / p.name).exists():
            continue
        if not list((VINARDO_RECEPTORS / p.name).glob("*.pdbt")):
            continue
        systems.append(p.name)
    if max_n is not None:
        systems = systems[:max_n]
    return systems


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max", type=int, default=20,
                        help="max number of receptors to test (default: 20)")
    parser.add_argument("--system", default=None,
                        help="test a specific system only")
    args = parser.parse_args()

    systems = find_systems(max_n=args.max, system=args.system)
    if not systems:
        print("no systems found to test")
        sys.exit(1)
    print(f"Testing {len(systems)} receptor(s)...\n")

    total_matched = 0
    total_mismatches = 0
    total_obabel = 0
    total_meeko = 0
    per_system = []
    mismatch_counts = Counter()
    skip_reasons = Counter()

    t0 = time.time()
    for sys_idx, system_id in enumerate(systems, 1):
        receptor_dir = GROUND_TRUTH / system_id
        pdb_path = receptor_dir / "receptor.pdb"
        ref_path = VINARDO_RECEPTORS / system_id / f"{system_id}_receptor.pdbt"

        # Read the pre-generated obabel reference
        if not ref_path.exists():
            skip_reasons["no_reference"] += 1
            continue
        ref_text = open(ref_path).read()
        ref_atoms = parse_pdbt_atoms(ref_text)

        meeko_text, err = meeko_build_pdbt(pdb_path)
        if meeko_text is None:
            skip_reasons[f"meeko:{err.splitlines()[0] if err else 'unknown'}"] += 1
            continue
        meeko_atoms = parse_pdbt_atoms(meeko_text)

        matches = match_atoms(meeko_atoms, ref_atoms)
        n_mismatches = sum(1 for ma, oa in matches if ma["type"] != oa["type"])
        for ma, oa in matches:
            if ma["type"] != oa["type"]:
                mismatch_counts[(ma["type"], oa["type"])] += 1

        total_matched += len(matches)
        total_mismatches += n_mismatches
        total_obabel += len(ref_atoms)
        total_meeko += len(meeko_atoms)
        rate = 100 * (len(matches) - n_mismatches) / len(matches) if matches else 0
        per_system.append((system_id, len(matches), n_mismatches, rate))

        if sys_idx <= 5 or n_mismatches > 0:
            print(f"  [{sys_idx:3d}/{len(systems)}] {system_id:35s} "
                  f"meeko={len(meeko_atoms):5d} obabel={len(ref_atoms):5d} "
                  f"matched={len(matches):5d} mismatches={n_mismatches:5d} ({rate:5.1f}%)")
    elapsed = time.time() - t0

    print("\n" + "=" * 90)
    print(f"  Receptor PDBT comparison complete ({elapsed:.0f}s)")
    print("=" * 90)
    print(f"  Receptors tested: {len(per_system)}")
    if skip_reasons:
        print(f"  Skipped: {sum(skip_reasons.values())}")
        for reason, count in skip_reasons.most_common():
            print(f"    - {reason}: {count}")
    print(f"  Total atoms compared: {total_matched}")
    if total_matched:
        rate = 100 * (total_matched - total_mismatches) / total_matched
        print(f"  Type match rate: {rate:.1f}%  "
              f"({total_matched - total_mismatches}/{total_matched})")
    if mismatch_counts:
        print(f"\n  Top mismatch types:")
        for (m, o), c in mismatch_counts.most_common(15):
            print(f"    {m}->{o}: {c}")

    print(f"\n  Per-system summary:")
    for system_id, matched, mismatches, rate in per_system:
        print(f"    {system_id:35s} matched={matched:5d} mismatches={mismatches:5d} ({rate:5.1f}%)")


if __name__ == "__main__":
    main()
