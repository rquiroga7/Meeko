"""Batch PDBT comparison: meeko-pdbt vs obabel-25-07 reference PDBTs.

Iterates over all ligands in vinardo_inputs/ligands/ with corresponding
SDF files in ground_truth/<system>/ligand_files/.
"""

import os
import sys
import math
import time
from pathlib import Path
from collections import Counter

from rdkit import Chem
from meeko import MoleculePreparation
from meeko.pdbt_writer import PDBTWriterLegacy
from meeko.pdbt_atomtyper import assign_pdbt_types
from rdkit.Chem import AllChem

REPO_ROOT = Path(__file__).resolve().parents[1]
DATASET = Path("/home/rquiroga/Datasets/runs-n-poses-datasets")
VINARDO_LIGANDS = DATASET / "vinardo_inputs" / "ligands"
GROUND_TRUTH = DATASET / "ground_truth"

USE_SYSTEM_OBABEL = True  # use system obabel-25-07 from PATH or custom path
OBABEL = REPO_ROOT / "openbabel-25-07" / "build-static" / "bin" / "obabel"

SUMMARY = {
    "total": 0,
    "errors": 0,
    "no_sdf": 0,
    "matched_atoms": 0,
    "total_mismatches": 0,
    "total_charge_diffs": 0,
    "mismatch_counts": Counter(),
}


def parse_pdbt_atoms(pdbt_text):
    atoms = []
    for line in pdbt_text.splitlines():
        if line.startswith("ATOM") or line.startswith("HETATM"):
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            charge_str = line[70:76].strip()
            atom_type = line[77:79].strip()
            try:
                charge = float(charge_str) if charge_str else 0.0
            except ValueError:
                charge = 0.0
            atoms.append({
                "coord": (x, y, z),
                "charge": charge,
                "type": atom_type,
            })
    return atoms


def _is_hydrogen_type(t):
    return t and t[0] == "H"

def match_atoms(meeko_atoms, obabel_atoms):
    """Two-phase matching: tight tolerance for heavy atoms, wider for hydrogens."""
    meeko_matched = [False] * len(meeko_atoms)
    obabel_matched = [False] * len(obabel_atoms)
    matches = []

    # Phase 1: heavy atoms (type starts with H -> hydrogen)
    for mi, ma in enumerate(meeko_atoms):
        if _is_hydrogen_type(ma["type"]):
            continue
        best_dist = 0.01
        best_oi = None
        for oi, oa in enumerate(obabel_atoms):
            if obabel_matched[oi]:
                continue
            if _is_hydrogen_type(oa["type"]):
                continue
            dist = math.sqrt(sum((ma["coord"][j] - oa["coord"][j])**2 for j in range(3)))
            if dist < best_dist:
                best_dist = dist
                best_oi = oi
        if best_oi is not None:
            meeko_matched[mi] = True
            obabel_matched[best_oi] = True
            matches.append((ma, obabel_atoms[best_oi]))

    # Phase 2: hydrogens with wider tolerance
    for mi, ma in enumerate(meeko_atoms):
        if meeko_matched[mi]:
            continue
        best_dist = 0.5
        best_oi = None
        for oi, oa in enumerate(obabel_atoms):
            if obabel_matched[oi]:
                continue
            dist = math.sqrt(sum((ma["coord"][j] - oa["coord"][j])**2 for j in range(3)))
            if dist < best_dist:
                best_dist = dist
                best_oi = oi
        if best_oi is not None:
            meeko_matched[mi] = True
            obabel_matched[best_oi] = True
            matches.append((ma, obabel_atoms[best_oi]))

    return matches


def meeko_prepare_ligand(sdf_text):
    mol = Chem.MolFromMolBlock(sdf_text, removeHs=False)
    if mol is None:
        return None
    mol = Chem.AddHs(mol, addCoords=True)
    try:
        preparator = MoleculePreparation()
        molsetups = preparator.prepare(mol)
    except Exception:
        return None
    if not molsetups:
        return None
    molsetup = molsetups[0]
    try:
        assign_pdbt_types(molsetup, mol)
    except Exception:
        return None
    pdbt_string, success, error_msg = PDBTWriterLegacy.write_string(molsetup)
    if not success:
        return None
    return pdbt_string


def compare_one(system_id):
    ligand_dir = GROUND_TRUTH / system_id / "ligand_files"
    if not ligand_dir.exists():
        return None, "no ligand_files dir"
    sdf_files = list(ligand_dir.glob("*.sdf"))
    if not sdf_files:
        return None, "no SDF files"
    sdf_path = sdf_files[0]

    ref_dir = VINARDO_LIGANDS / system_id
    if not ref_dir.exists():
        return None, "no reference dir"
    pdbt_files = list(ref_dir.glob("*.pdbt"))
    if not pdbt_files:
        return None, "no reference PDBT"

    # Read reference PDBT
    ref_text = open(pdbt_files[0]).read()
    ref_atoms = parse_pdbt_atoms(ref_text)

    # Run meeko
    sdf_text = open(sdf_path).read()
    meeko_text = meeko_prepare_ligand(sdf_text)
    if meeko_text is None:
        return None, "meeko failed"
    meeko_atoms = parse_pdbt_atoms(meeko_text)

    # Match
    matches = match_atoms(meeko_atoms, ref_atoms)

    type_mismatches = []
    charge_diffs = 0
    for ma, oa in matches:
        if ma["type"] != oa["type"]:
            type_mismatches.append((ma["type"], oa["type"]))
        if abs(ma["charge"] - oa["charge"]) > 0.001:
            charge_diffs += 1

    return {
        "system": system_id,
        "meeko_atoms": len(meeko_atoms),
        "ref_atoms": len(ref_atoms),
        "matched": len(matches),
        "type_mismatches": type_mismatches,
        "charge_diffs": charge_diffs,
    }, None


def format_time(sec):
    if sec < 60:
        return f"{sec:.0f}s"
    return f"{sec//60}m{sec%60:02d}s"


def main():
    # Collect all systems that have both ground truth and vinardo references
    systems = []
    for d in sorted(VINARDO_LIGANDS.iterdir()):
        if d.is_dir() and (GROUND_TRUTH / d.name / "ligand_files").exists():
            systems.append(d.name)

    total = len(systems)
    print(f"Found {total} systems with both SDF and reference PDBT")
    print(f"Meeko: {REPO_ROOT}")
    print(f"Dataset: {DATASET}")
    print(f"Obabel ref: {VINARDO_LIGANDS}")
    print()
    print(f"{'Progress':>8s}  {'System':>40s}  {'Matched':>7s}  {'OK':>5s}  {'Diff':>5s}  {'ChargeΔ':>7s}  {'Time':>8s}")
    print("-" * 90)

    start_time = time.time()
    t0 = start_time

    ok_atoms = 0
    diff_atoms = 0
    charge_diff_count = 0
    error_count = 0

    for i, system_id in enumerate(systems):
        result, err = compare_one(system_id)
        t_now = time.time()
        elapsed = t_now - t0

        if result is None:
            error_count += 1
            print(f"  {i+1:5d}/{total}  {system_id:40s}  {'ERROR':>7s}  {'':>5s}  {'':>5s}  {err:>20s}  {format_time(elapsed):>8s}")
            continue

        n_mismatch = len(result["type_mismatches"])
        n_charge = result["charge_diffs"]
        matched = result["matched"]
        n_ok = matched - n_mismatch

        ok_atoms += n_ok
        diff_atoms += n_mismatch
        charge_diff_count += n_charge
        SUMMARY["total_mismatches"] += n_mismatch
        for mt, mo in result["type_mismatches"]:
            SUMMARY["mismatch_counts"][f"{mt}->{mo}"] += 1

        status = "OK" if n_mismatch == 0 else f"  {n_mismatch:2d}"
        print(f"  {i+1:5d}/{total}  {system_id:40s}  {matched:7d}  {n_ok:5d}  {n_mismatch:5d}  {n_charge:7d}  {format_time(elapsed):>8s}")

        SUMMARY["matched_atoms"] += matched

    total_time = time.time() - start_time

    print("-" * 90)
    print(f"\n{'='*70}")
    print(f"  PDBT Batch Validation Complete")
    print(f"  Systems: {total}  Errors: {error_count}")
    print(f"  Total atoms matched: {ok_atoms + diff_atoms}")
    print(f"  Type matches: {ok_atoms}  Type mismatches: {diff_atoms}")
    print(f"  Overall type match rate: {ok_atoms / max(ok_atoms + diff_atoms, 1) * 100:.1f}%")
    print(f"  Ligands with at least 1 charge diff: {charge_diff_count}")
    print(f"  Total time: {format_time(total_time)}")
    print(f"{'='*70}")

    if diff_atoms > 0:
        print(f"\n  Top mismatch types:")
        for pair, count in SUMMARY["mismatch_counts"].most_common(30):
            print(f"    {pair:>10s}: {count:5d}")


if __name__ == "__main__":
    main()
