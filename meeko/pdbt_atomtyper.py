from rdkit import Chem

_AMINOACID_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
    "HIS", "HID", "HIE", "HIP", "ILE", "LEU", "LYS", "MET",
    "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "XAA", "ASX", "XLE", "GLX", "CYX", "SEC", "PYL",
}

# Hardcoded PDBT atom types for standard protein residues, matching OpenBabel-25-07.
# Backbone atoms (C, N, O, OXT) apply to all residues; CA uses GLY vs non-GLY rules.
# Atoms not found here fall through to chemical-environment typing.
_PROTEIN_BACKBONE = {
    "C": "A1",
    "N": "Nf",
    "O": "Of",
    "OXT": "Oc",
}

# Sidechain atoms per residue, extracted from openbabel-25-07/src/formats/pdbtformat.cpp:1408-1586
# HIS variants:
#   HIS  → ND1=Np (acceptor), NE2=Nb (donor)
#   HID  → ND1=Nb (donor),    NE2=Np (acceptor)
#   HIE  → ND1=Nb (donor),    NE2=Nb (donor)
#   HIP  → ND1=Nb (donor),    NE2=Nb (donor)
_PROTEIN_ATOM_TYPES = {
    "ALA": {"CB": "C3"},
    "ARG": {"CB": "C2", "CG": "C2", "CD": "C2", "NE": "Nb", "CZ": "C1", "NH1": "Na", "NH2": "Na"},
    "ASN": {"CB": "C2", "CG": "A1", "OD1": "Of", "ND2": "Nf"},
    "ASP": {"CB": "C2", "CG": "A1", "OD1": "Oc", "OD2": "Oc"},
    "CYS": {"CB": "C2", "SG": "S1"},
    "GLN": {"CB": "C2", "CG": "C2", "CD": "A1", "OE1": "Of", "NE2": "Nf"},
    "GLU": {"CB": "C2", "CG": "C2", "CD": "A1", "OE1": "Oc", "OE2": "Oc"},
    "HIS": {"CB": "C2", "CG": "A0", "CD2": "A0", "CE1": "A0", "ND1": "Np", "NE2": "Nb"},
    "HID": {"CB": "C2", "CG": "A0", "CD2": "A0", "CE1": "A0", "ND1": "Nb", "NE2": "Np"},
    "HIE": {"CB": "C2", "CG": "A0", "CD2": "A0", "CE1": "A0", "ND1": "Nb", "NE2": "Nb"},
    "HIP": {"CB": "C2", "CG": "A0", "CD2": "A0", "CE1": "A0", "ND1": "Nb", "NE2": "Nb"},
    "ILE": {"CB": "C1", "CG1": "C2", "CG2": "C3", "CD1": "C3"},
    "LEU": {"CB": "C2", "CG": "C1", "CD1": "C3", "CD2": "C3"},
    "LYS": {"CB": "C2", "CG": "C2", "CD": "C2", "CE": "C2", "NZ": "Ni"},
    "MET": {"CB": "C2", "CG": "C2", "SD": "S3", "CE": "C3"},
    "PHE": {"CB": "C2", "CG": "A0", "CD1": "A0", "CD2": "A0", "CE1": "A0", "CE2": "A0", "CZ": "A0"},
    "PRO": {"CB": "C2", "CG": "C2", "CD": "C2"},
    "SER": {"CB": "C2", "OG": "Ob"},
    "THR": {"CB": "C1", "OG1": "Ob", "CG2": "C3"},
    "TRP": {"CB": "C2", "CG": "A0", "CD1": "A0", "CD2": "A0", "CE2": "A0", "CE3": "A0", "CZ2": "A0", "CZ3": "A0", "CH2": "A0", "NE1": "Nb"},
    "TYR": {"CB": "C2", "CG": "A0", "CD1": "A0", "CD2": "A0", "CE1": "A0", "CE2": "A0", "CZ": "A0", "OH": "Oa"},
    "VAL": {"CB": "C1", "CG1": "C3", "CG2": "C3"},
}


def _lookup_protein_pdbt_type(res_name, atom_name):
    """Look up hardcoded PDBT atom type for a standard amino-acid atom.

    Returns None if no hardcoded type exists (fall through to chemical rules).
    """
    res_upper = res_name.upper()[:3]
    atom_clean = atom_name.strip()

    if atom_clean in _PROTEIN_BACKBONE:
        return _PROTEIN_BACKBONE[atom_clean]

    if atom_clean == "CA":
        return "C2" if res_upper == "GLY" else "C1"

    res_table = _PROTEIN_ATOM_TYPES.get(res_upper)
    if res_table is not None:
        return res_table.get(atom_clean)

    return None


def _hvy_degree(atom):
    return sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() != 1)


def _count_terminal(atom, atomic_num):
    count = 0
    for nbr in atom.GetNeighbors():
        if nbr.GetAtomicNum() == atomic_num and _hvy_degree(nbr) == 1:
            count += 1
    return count


def _count_free_oxygens(atom):
    return _count_terminal(atom, 8)


def _count_free_nitrogens(atom):
    return _count_terminal(atom, 7)


# ---------------------------------------------------------------------------
# Ring / aromaticity helpers
# ---------------------------------------------------------------------------

def _is_in_ring_poor_def(atom):
    if atom.GetIsAromatic() and atom.IsInRing():
        return True
    if atom.GetDegree() < 2:
        return False
    if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP2:
        return False
    if not atom.IsInRing():
        return False
    count = 0
    for nbr in atom.GetNeighbors():
        if nbr.GetAtomicNum() == 1:
            continue
        if nbr.GetHybridization() == Chem.rdchem.HybridizationType.SP2 and nbr.IsInRing():
            count += 1
    return count >= 2


# ---------------------------------------------------------------------------
# H-bond helpers
# ---------------------------------------------------------------------------

def _is_hbond_donor(atom):
    if atom.GetAtomicNum() not in (7, 8, 9, 16):
        return False
    for nbr in atom.GetNeighbors():
        if nbr.GetAtomicNum() == 1:
            return True
    return False


def _is_hbond_acceptor(atom):
    num = atom.GetAtomicNum()
    if num == 7:
        if atom.GetFormalCharge() > 0:
            return False
        degree = atom.GetDegree()
        if degree < 4:
            return True
        # N with degree 4 (e.g. ammonium) is not an acceptor
        return False
    elif num == 8:
        if atom.GetDegree() < 3:
            return True
        return False
    elif num == 16:
        if atom.GetDegree() < 3:
            return True
        return False
    return False


# ---------------------------------------------------------------------------
# Carbon helpers
# ---------------------------------------------------------------------------

def _is_carbonyl_oxygen(atom, mol):
    if atom.GetAtomicNum() != 8:
        return False
    if _hvy_degree(atom) != 1:
        return False
    for nbr in atom.GetNeighbors():
        if nbr.GetAtomicNum() != 6:
            continue
        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
        if bond is not None and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            return True
    return False


def _is_amide_bond(bond, mol):
    begin = bond.GetBeginAtom()
    end = bond.GetEndAtom()
    if begin.GetAtomicNum() == 6 and end.GetAtomicNum() == 7:
        c_atom, n_atom = begin, end
    elif begin.GetAtomicNum() == 7 and end.GetAtomicNum() == 6:
        c_atom, n_atom = end, begin
    else:
        return False
    if c_atom.GetHybridization() != Chem.rdchem.HybridizationType.SP2:
        return False
    for nbr in c_atom.GetNeighbors():
        nbr_bond = mol.GetBondBetweenAtoms(c_atom.GetIdx(), nbr.GetIdx())
        if nbr_bond is not None and nbr.GetAtomicNum() == 8 and nbr_bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            return True
    return False


def _is_carboxyl_oxygen(atom):
    if atom.GetAtomicNum() != 8:
        return False
    if _hvy_degree(atom) != 1:
        return False
    for nbr in atom.GetNeighbors():
        if nbr.GetAtomicNum() == 6 and _count_terminal(nbr, 8) >= 2:
            return True
    return False


def _is_phosphate_oxygen(atom):
    if atom.GetAtomicNum() != 8:
        return False
    if _hvy_degree(atom) != 1:
        return False
    for nbr in atom.GetNeighbors():
        if nbr.GetAtomicNum() == 15 and _count_free_oxygens(nbr) > 1:
            return True
    return False


def _is_phosphate_bridge_oxygen(atom):
    if atom.GetAtomicNum() != 8:
        return False
    if _hvy_degree(atom) != 2:
        return False
    for nbr in atom.GetNeighbors():
        if nbr.GetAtomicNum() == 15 and _count_free_oxygens(nbr) > 1:
            return True
    return False


def _is_ether_oxygen(atom):
    if atom.GetAtomicNum() != 8:
        return False
    if _hvy_degree(atom) != 2:
        return False
    carbon_count = sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() == 6)
    return carbon_count == 2


def _is_phenol_oxygen(atom):
    if atom.GetAtomicNum() != 8:
        return False
    if _hvy_degree(atom) != 1:
        return False
    for nbr in atom.GetNeighbors():
        if nbr.GetAtomicNum() == 6 and nbr.GetIsAromatic():
            return True
    return False


def _is_sulfate_oxygen(atom):
    if atom.GetAtomicNum() != 8:
        return False
    if _hvy_degree(atom) != 1:
        return False
    for nbr in atom.GetNeighbors():
        if nbr.GetAtomicNum() == 16 and _count_terminal(nbr, 8) >= 3:
            return True
    return False


def _is_nitro_oxygen(atom):
    if atom.GetAtomicNum() != 8:
        return False
    if _hvy_degree(atom) != 1:
        return False
    for nbr in atom.GetNeighbors():
        if nbr.GetAtomicNum() == 7 and _count_terminal(nbr, 8) >= 2:
            return True
    return False


# ---------------------------------------------------------------------------
# Amide / ester helpers (O-oriented)
# ---------------------------------------------------------------------------

def _is_amide_oxygen(atom, mol, aromatic):
    if atom.GetAtomicNum() != 8:
        return False
    if _hvy_degree(atom) != 1:
        return False
    for nbr in atom.GetNeighbors():
        if nbr.GetAtomicNum() == 6:
            if nbr.GetIsAromatic() != aromatic:
                continue
            if nbr.GetHybridization() != Chem.rdchem.HybridizationType.SP2:
                continue
            for bond in mol.GetBonds():
                if bond.GetBeginAtomIdx() == nbr.GetIdx() or bond.GetEndAtomIdx() == nbr.GetIdx():
                    if _is_amide_bond(bond, mol):
                        return True
    return False


def _is_amide_aroma_oxygen(atom, mol):
    return _is_amide_oxygen(atom, mol, True)


def _is_amide_alifa_oxygen(atom, mol):
    return _is_amide_oxygen(atom, mol, False)


def _is_ester_carbonyl_oxygen(atom, mol):
    if atom.GetAtomicNum() != 8:
        return False
    if _hvy_degree(atom) != 1:
        return False
    for nbr in atom.GetNeighbors():
        if nbr.GetAtomicNum() == 6:
            for bond in mol.GetBonds():
                if bond.GetBeginAtomIdx() == nbr.GetIdx() or bond.GetEndAtomIdx() == nbr.GetIdx():
                    if _is_ester_bond(bond, mol):
                        return True
    return False


def _is_ester_oxygen(atom, mol):
    if atom.GetAtomicNum() != 8:
        return False
    if _hvy_degree(atom) != 2:
        return False
    for bond in mol.GetBonds():
        if bond.GetBeginAtomIdx() == atom.GetIdx() or bond.GetEndAtomIdx() == atom.GetIdx():
            if _is_ester_bond(bond, mol):
                return True
    return False


def _is_ester_bond(bond, mol):
    begin = bond.GetBeginAtom()
    end = bond.GetEndAtom()
    if begin.GetAtomicNum() == 6 and end.GetAtomicNum() == 8:
        c_atom = begin
    elif begin.GetAtomicNum() == 8 and end.GetAtomicNum() == 6:
        c_atom = end
    else:
        return False
    if c_atom.GetHybridization() != Chem.rdchem.HybridizationType.SP2:
        return False
    for nbr in c_atom.GetNeighbors():
        nbr_bond = mol.GetBondBetweenAtoms(c_atom.GetIdx(), nbr.GetIdx())
        if nbr_bond is not None and nbr.GetAtomicNum() == 8 and nbr_bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            return True
    return False


# ---------------------------------------------------------------------------
# Sulfur helpers
# ---------------------------------------------------------------------------

def _is_sulfone_oxygen(atom):
    if atom.GetAtomicNum() != 8:
        return False
    if _hvy_degree(atom) != 1:
        return False
    for nbr in atom.GetNeighbors():
        if nbr.GetAtomicNum() == 16 and _count_free_oxygens(nbr) == 2 and _hvy_degree(nbr) == 4:
            carbon_count = sum(1 for n2 in nbr.GetNeighbors() if n2.GetAtomicNum() == 6)
            if carbon_count == 2:
                return True
    return False


def _is_sulfonamide_oxygen(atom):
    if atom.GetAtomicNum() != 8:
        return False
    if _hvy_degree(atom) != 1:
        return False
    for nbr in atom.GetNeighbors():
        if nbr.GetAtomicNum() == 16 and _count_free_oxygens(nbr) == 2:
            for n2 in nbr.GetNeighbors():
                if n2.GetAtomicNum() == 7:
                    return True
    return False


def _is_sulfonamide_nitrogen(atom):
    if atom.GetAtomicNum() != 7:
        return False
    for nbr in atom.GetNeighbors():
        if nbr.GetAtomicNum() == 16:
            for n2 in nbr.GetNeighbors():
                if n2.GetAtomicNum() == 8 and _hvy_degree(n2) == 1:
                    return True
    return False


def _is_sulfur_acceptor(atom):
    return atom.GetAtomicNum() == 16 and _hvy_degree(atom) <= 2


# ---------------------------------------------------------------------------
# Nitrogen helpers
# ---------------------------------------------------------------------------

def _is_aniline_nitrogen(atom):
    if atom.GetAtomicNum() != 7:
        return False
    if _hvy_degree(atom) != 1:
        return False
    count_h = 0
    count_a = 0
    for nbr in atom.GetNeighbors():
        if nbr.GetAtomicNum() == 1:
            count_h += 1
        if nbr.GetAtomicNum() == 6 and _is_in_ring_poor_def(nbr):
            count_a += 1
    return count_h == 2 and count_a == 1


def _is_aniline_sus_nitrogen(atom):
    if atom.GetAtomicNum() != 7:
        return False
    if _hvy_degree(atom) != 2:
        return False
    count_h = 0
    count_c = 0
    count_a = 0
    for nbr in atom.GetNeighbors():
        if nbr.GetAtomicNum() == 1:
            count_h += 1
        if nbr.GetAtomicNum() == 6 and _is_in_ring_poor_def(nbr):
            count_a += 1
        if nbr.GetAtomicNum() == 6 and nbr.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
            count_c += 1
    return count_h == 1 and count_c == 1 and count_a == 1


def _is_guanidinium_nitrogen(atom):
    if atom.GetAtomicNum() != 7:
        return False
    for nbr in atom.GetNeighbors():
        if nbr.GetAtomicNum() == 6 and _hvy_degree(nbr) == 3 and _count_free_nitrogens(nbr) == 2:
            if sum(1 for n2 in nbr.GetNeighbors() if n2.GetAtomicNum() == 7) == 3:
                return True
    return False


def _is_amidine_nitrogen(atom):
    if atom.GetAtomicNum() != 7:
        return False
    if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP2:
        return False
    if _hvy_degree(atom) != 1:
        return False
    for nbr in atom.GetNeighbors():
        if nbr.GetAtomicNum() == 6 and _hvy_degree(nbr) == 3 and _count_free_nitrogens(nbr) == 2:
            nitrogen_count = 0
            carbon_count = 0
            for n2 in nbr.GetNeighbors():
                if n2.GetAtomicNum() == 7:
                    nitrogen_count += 1
                if n2.GetAtomicNum() == 6:
                    carbon_count += 1
            if nitrogen_count == 2 and carbon_count == 1:
                return True
    return False


def _is_nitrile_nitrogen(atom):
    if atom.GetAtomicNum() != 7:
        return False
    if _hvy_degree(atom) != 1 or atom.GetHybridization() != Chem.rdchem.HybridizationType.SP:
        return False
    for nbr in atom.GetNeighbors():
        if nbr.GetAtomicNum() == 6 and nbr.GetHybridization() == Chem.rdchem.HybridizationType.SP and _hvy_degree(nbr) == 2:
            return True
    return False


def _is_nitrogen_charged_sp3(atom):
    if atom.GetAtomicNum() != 7:
        return False
    if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
        return False
    carbon_sp3_count = 0
    for nbr in atom.GetNeighbors():
        if nbr.GetAtomicNum() == 6 and nbr.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
            carbon_sp3_count += 1
        elif nbr.GetAtomicNum() == 1:
            continue
        else:
            return False
    return 0 < carbon_sp3_count < 4


def _is_aromatic_nitrogen(atom):
    if atom.GetAtomicNum() != 7:
        return False
    if _is_in_ring_poor_def(atom):
        return True
    count = 0
    for nbr in atom.GetNeighbors():
        if nbr.GetAtomicNum() == 1:
            continue
        if _is_in_ring_poor_def(nbr):
            count += 1
    return count == 2


def _is_amide_nitrogen(atom, mol):
    if atom.GetAtomicNum() != 7:
        return False
    for nbr in atom.GetNeighbors():
        if nbr.GetAtomicNum() == 6:
            for bond in mol.GetBonds():
                if (bond.GetBeginAtomIdx() == nbr.GetIdx() or bond.GetEndAtomIdx() == nbr.GetIdx()):
                    if _is_amide_bond(bond, mol):
                        return True
    return False


def _is_aminoacid_residue(res_name):
    return res_name.upper()[:3] in _AMINOACID_RESIDUES


# ---------------------------------------------------------------------------
# Main typing functions
# ---------------------------------------------------------------------------

def _type_carbon(atom, mol):
    if atom.GetIsAromatic():
        return "A0"
    hyb = atom.GetHybridization()
    if hyb == Chem.rdchem.HybridizationType.SP3:
        hd = _hvy_degree(atom)
        if hd == 0:
            return "C1"
        elif hd == 1:
            return "C3"
        elif hd == 2:
            return "C2"
        else:
            return "C1"
    elif hyb == Chem.rdchem.HybridizationType.SP2:
        if atom.IsInRing():
            return "A0"
        else:
            return "A1"
    elif hyb == Chem.rdchem.HybridizationType.SP:
        return "A1"
    else:
        return "C1"


def _type_oxygen(atom, mol):
    if _is_hbond_donor(atom) and _is_hbond_acceptor(atom):
        if _is_phosphate_oxygen(atom):
            return "Od"
        if _is_carboxyl_oxygen(atom):
            return "Oc"
        if _is_sulfate_oxygen(atom):
            return "Oj"
        if _is_nitro_oxygen(atom):
            return "Ok"
        if _is_sulfonamide_oxygen(atom):
            return "Ol"
        if _is_phenol_oxygen(atom):
            return "Oa"
        else:
            return "Ob"
    if _hvy_degree(atom) == 1:
        if _is_carboxyl_oxygen(atom):
            return "Oc"
        if _is_phosphate_oxygen(atom):
            return "Od"
        if _is_amide_aroma_oxygen(atom, mol):
            return "Oe"
        if _is_amide_alifa_oxygen(atom, mol):
            return "Of"
        if _is_ester_carbonyl_oxygen(atom, mol):
            return "Og"
        if _is_carbonyl_oxygen(atom, mol):
            return "Oh"
        if _is_sulfone_oxygen(atom):
            return "Oi"
        if _is_sulfate_oxygen(atom):
            return "Oj"
        if _is_nitro_oxygen(atom):
            return "Ok"
        if _is_sulfonamide_oxygen(atom):
            return "Ol"
        return "Om"
    elif _hvy_degree(atom) == 2:
        if _is_phosphate_bridge_oxygen(atom):
            return "On"
        if _is_ester_oxygen(atom, mol):
            return "Oo"
        if _is_ether_oxygen(atom):
            return "Op"
        if atom.GetIsAromatic():
            return "Oq"
        return "Or"
    return "Ox"


def _type_nitrogen(atom, mol):
    if _is_hbond_donor(atom):
        if atom.GetHybridization() == Chem.rdchem.HybridizationType.SP2:
            if _is_guanidinium_nitrogen(atom):
                return "Na"
            if _is_aromatic_nitrogen(atom):
                return "Nb"
            if _is_sulfonamide_nitrogen(atom):
                return "Nc"
            if _is_aniline_nitrogen(atom):
                return "Nd"
            if _is_aniline_sus_nitrogen(atom):
                return "Ne"
            if _is_amide_nitrogen(atom, mol):
                return "Nf"
            if _is_amidine_nitrogen(atom):
                return "Ng"
            return "Nh"
        elif atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
            if _is_nitrogen_charged_sp3(atom):
                return "Ni"
            if _is_sulfonamide_nitrogen(atom):
                return "Nj"
            return "Nk"
        else:
            return "Nl"
    elif _is_hbond_acceptor(atom):
        if atom.GetHybridization() == Chem.rdchem.HybridizationType.SP:
            if _is_nitrile_nitrogen(atom):
                return "Nm"
            return "Nn"
        elif atom.GetHybridization() == Chem.rdchem.HybridizationType.SP2:
            if _is_amidine_nitrogen(atom):
                return "No"
            if _is_in_ring_poor_def(atom):
                hetero_degree = sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() not in (1, 6))
                if hetero_degree == 0:
                    return "Np"
                else:
                    return "Nq"
            return "Nr"
        elif atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
            if _is_nitrogen_charged_sp3(atom):
                return "Ns"
            return "Nt"
        else:
            return "Nt"
    else:
        if atom.GetHybridization() == Chem.rdchem.HybridizationType.SP2:
            return "Nu"
        elif atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
            return "Nv"
        else:
            return "Nu"


def _type_sulfur(atom, mol):
    if _is_hbond_donor(atom):
        return "S1"
    if _is_sulfur_acceptor(atom):
        return "S2"
    return "S3"


def _get_default_type(atom):
    symbol = atom.GetSymbol()
    if len(symbol) >= 2:
        return symbol[0].upper() + symbol[1].lower()
    return symbol + " "


def get_pdbt_atom_type(atom, mol):
    # Fast hardcoded lookup for known protein residues (avoids chemical analysis)
    pdb_info = atom.GetMonomerInfo()
    if pdb_info is not None and hasattr(pdb_info, "GetResidueName"):
        res_name = pdb_info.GetResidueName()
        atom_name = pdb_info.GetName()
        if _is_aminoacid_residue(res_name):
            pdbt_type = _lookup_protein_pdbt_type(res_name, atom_name)
            if pdbt_type is not None:
                return pdbt_type
    atomic_num = atom.GetAtomicNum()
    if atomic_num == 1:
        return "HD"
    if atomic_num == 6:
        return _type_carbon(atom, mol)
    if atomic_num == 8:
        return _type_oxygen(atom, mol)
    if atomic_num == 7:
        return _type_nitrogen(atom, mol)
    if atomic_num == 16:
        return _type_sulfur(atom, mol)
    return _get_default_type(atom)


def assign_pdbt_types(molsetup, mol):
    for atom_idx in range(molsetup.true_atom_count):
        atom = mol.GetAtomWithIdx(atom_idx)
        pdbt_type = get_pdbt_atom_type(atom, mol)
        molsetup.set_atom_type(atom_idx, pdbt_type)
