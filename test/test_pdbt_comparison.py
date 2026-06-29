"""Compare meeko-pdbt with obabel-25-07 PDBT output on runs-n-poses dataset.

PDBT format (1-indexed columns):
  Col  1- 6: record type (ATOM  / HETATM)
  Col  7-11: atom serial
  Col 12   : space
  Col 13-16: atom name (4 chars)
  Col 17   : altLoc
  Col 18-20: residue name (3 chars)
  Col 21   : space
  Col 22   : chain
  Col 23-26: residue number
  Col 27   : insertion code
  Col 28-30: spaces
  Col 31-38: x (8 chars)
  Col 39-46: y (8 chars)
  Col 47-54: z (8 chars)
  Col 55-59: occupancy / vdW (5 chars)
  Col 60-65: tempFactor / Elec (6 chars)
  Col 66-69: spaces (4 chars)
  Col 70-75: charge (6 chars, e.g. +0.000)
  Col 76   : space
  Col 77-78: atom type (2 chars)
"""

import os
import sys
import subprocess as sp
import tempfile
import math
from pathlib import Path

from rdkit import Chem
from meeko import MoleculePreparation
from meeko.pdbt_writer import PDBTWriterLegacy
from meeko.pdbt_atomtyper import assign_pdbt_types

REPO_ROOT = Path(__file__).resolve().parents[1]
RUNS_N_POSES = REPO_ROOT / "runs-n-poses"
DATASET = Path("/home/rquiroga/Datasets/runs-n-poses-datasets")
OBABEL = REPO_ROOT / "openbabel-25-07" / "build-static" / "bin" / "obabel"

TEST_SYSTEMS = [
    "5s9z__1__1.A_1.B__1.R",
    "5sau__1__1.A__1.B",
    "5s9y__1__1.A__1.K",
    "7nag__1__1.A__1.C",
]

PASS = 0
FAIL = 0
RESULTS = []


def parse_pdbt_atoms(pdbt_text):
    atoms = []
    for line in pdbt_text.splitlines():
        if line.startswith("ATOM") or line.startswith("HETATM"):
            serial = int(line[6:11].strip())
            name = line[12:16].strip()
            resname = line[17:20].strip()
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
                "serial": serial,
                "name": name,
                "resname": resname,
                "coord": (x, y, z),
                "charge": charge,
                "type": atom_type,
            })
    return atoms


def match_atoms(meeko_atoms, obabel_atoms, tolerance=0.01):
    """Match atoms between meeko and obabel by coordinate proximity."""
    meeko_matched = [False] * len(meeko_atoms)
    obabel_matched = [False] * len(obabel_atoms)
    matches = []

    for mi, ma in enumerate(meeko_atoms):
        best_dist = tolerance
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
            matches.append((ma, obabel_atoms[best_oi], best_dist))

    unmatched_meeko = [ma for mi, ma in enumerate(meeko_atoms) if not meeko_matched[mi]]
    unmatched_obabel = [oa for oi, oa in enumerate(obabel_atoms) if not obabel_matched[oi]]

    return matches, unmatched_meeko, unmatched_obabel


def compare_pdbt(meeko_text, obabel_text, label):
    meeko_atoms = parse_pdbt_atoms(meeko_text)
    obabel_atoms = parse_pdbt_atoms(obabel_text)

    matches, unmatched_meeko, unmatched_obabel = match_atoms(meeko_atoms, obabel_atoms)

    type_mismatches = []
    charge_diffs = []
    for ma, oa, dist in matches:
        if ma["type"] != oa["type"]:
            type_mismatches.append({
                "name": ma["name"],
                "coord": ma["coord"],
                "meeko_type": ma["type"],
                "obabel_type": oa["type"],
                "meeko_charge": ma["charge"],
                "obabel_charge": oa["charge"],
                "distance": dist,
            })
        if abs(ma["charge"] - oa["charge"]) > 0.001:
            charge_diffs.append({
                "name": ma["name"],
                "meeko_charge": ma["charge"],
                "obabel_charge": oa["charge"],
            })

    total_matched = len(matches)
    match_rate = (total_matched - len(type_mismatches)) / max(total_matched, 1) * 100

    return {
        "label": label,
        "meeko_atoms": len(meeko_atoms),
        "obabel_atoms": len(obabel_atoms),
        "matched": total_matched,
        "unmatched_meeko": unmatched_meeko,
        "unmatched_obabel": unmatched_obabel,
        "type_mismatches": type_mismatches,
        "charge_diffs": charge_diffs,
        "match_rate": match_rate,
    }


def meeko_prepare_ligand(sdf_path):
    mol = next(Chem.SDMolSupplier(str(sdf_path), removeHs=False))
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    preparator = MoleculePreparation()
    molsetups = preparator.prepare(mol)
    if not molsetups:
        return None
    molsetup = molsetups[0]
    assign_pdbt_types(molsetup, mol)
    pdbt_string, success, error_msg = PDBTWriterLegacy.write_string(molsetup)
    if not success:
        return None
    return pdbt_string


def test_ligand(system_id):
    ligand_dir = RUNS_N_POSES / system_id / "ligand_files"
    if not ligand_dir.exists():
        return None
    sdf_files = list(ligand_dir.glob("*.sdf"))
    if not sdf_files:
        return None
    sdf_path = sdf_files[0]
    chain = sdf_path.stem

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_sdf = os.path.join(tmpdir, f"{chain}.sdf")
        with open(tmp_sdf, "w") as f:
            f.write(open(sdf_path).read())

        obabel_pdbt = os.path.join(tmpdir, f"{chain}_obabel.pdbt")
        cmd = [str(OBABEL), "-i", "sdf", tmp_sdf, "-o", "pdbt", "-p7", "-O", obabel_pdbt]
        r = sp.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            return {"error": f"obabel failed: {r.stderr}"}

        with open(obabel_pdbt) as f:
            obabel_text = f.read()

    meeko_text = meeko_prepare_ligand(sdf_path)
    if meeko_text is None:
        return {"error": "meeko preparation failed"}

    return compare_pdbt(meeko_text, obabel_text, f"ligand:{system_id}/{chain}")


def run_tests():
    global PASS, FAIL
    print("=" * 70)
    print("PDBT Comparison: meeko-pdbt vs obabel-25-07 (runs-n-poses)")
    print("Coordinate-based atom matching")
    print("=" * 70)

    if not RUNS_N_POSES.exists():
        print(f"  SKIP  all: runs-n-poses symlink not found at {RUNS_N_POSES}")
        return False

    for system_id in TEST_SYSTEMS:
        result = test_ligand(system_id)
        if result is None:
            print(f"  SKIP  {system_id}: no ligand SDF found")
            continue
        if "error" in result:
            print(f"  FAIL  {system_id}: {result['error']}")
            FAIL += 1
            continue

        label = result["label"]
        n_mismatch = len(result["type_mismatches"])
        n_charge = len(result["charge_diffs"])
        n_meeko_un = len(result["unmatched_meeko"])
        n_obabel_un = len(result["unmatched_obabel"])

        print(f"\n  {'='*58}")
        print(f"  {label}")
        print(f"  {'='*58}")
        print(f"    Atoms: meeko={result['meeko_atoms']}, "
              f"obabel={result['obabel_atoms']}, matched={result['matched']}")
        print(f"    Unmatched meeko: {n_meeko_un}, obabel: {n_obabel_un}")
        print(f"    Type match rate: {result['match_rate']:.1f}%")
        print(f"    Type mismatches: {n_mismatch}")
        print(f"    Charge diffs (|dq|>0.001): {n_charge}")

        if n_meeko_un > 0:
            print(f"    Unmatched meeko atoms (up to 5):")
            for a in result["unmatched_meeko"][:5]:
                print(f"      {a['name']:>4s} @ ({a['coord'][0]:.2f}, {a['coord'][1]:.2f}, {a['coord'][2]:.2f})")
        if n_obabel_un > 0:
            print(f"    Unmatched obabel atoms (up to 5):")
            for a in result["unmatched_obabel"][:5]:
                print(f"      {a['name']:>4s} @ ({a['coord'][0]:.2f}, {a['coord'][1]:.2f}, {a['coord'][2]:.2f})")

        if n_mismatch > 0:
            print(f"    Type mismatches (up to 15):")
            for m in result["type_mismatches"][:15]:
                print(f"      {m['name']:>4s} @ ({m['coord'][0]:.2f}, {m['coord'][1]:.2f}, {m['coord'][2]:.2f}):  "
                      f"meeko={m['meeko_type']:>4s}  obabel={m['obabel_type']:>4s}  "
                      f"q_m={m['meeko_charge']:+.3f}  q_o={m['obabel_charge']:+.3f}")
            if len(result["type_mismatches"]) > 15:
                print(f"      ... and {len(result['type_mismatches']) - 15} more")

        PASS += 1
        RESULTS.append(result)

    print()
    print("=" * 70)
    total = PASS + FAIL
    print(f"Results: {PASS}/{total} passed, {FAIL}/{total} failed")
    print("=" * 70)

    if RESULTS:
        all_matched = sum(r["matched"] for r in RESULTS)
        all_mismatches = sum(len(r["type_mismatches"]) for r in RESULTS)
        all_charge = sum(len(r["charge_diffs"]) for r in RESULTS)
        print(f"Total atoms matched: {all_matched}")
        print(f"Total type mismatches: {all_mismatches}")
        print(f"Total charge diffs: {all_charge}")
        if all_matched > 0:
            print(f"Overall type match rate: "
                  f"{(all_matched - all_mismatches) / all_matched * 100:.1f}%")
    return FAIL == 0


if __name__ == "__main__":
    success = run_tests()
    sys.exit(0 if success else 1)
