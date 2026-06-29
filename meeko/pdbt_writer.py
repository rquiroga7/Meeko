import logging
import math
from rdkit import Chem
from .utils import pdbutils

logger = logging.getLogger(__name__)

_ELEMENT_SYMBOLS = {
    1: "H", 5: "B", 6: "C", 7: "N", 8: "O", 9: "F",
    12: "Mg", 14: "Si", 15: "P", 16: "S", 17: "Cl",
    20: "Ca", 25: "Mn", 26: "Fe", 30: "Zn", 35: "Br", 53: "I",
}


class PDBTWriterLegacy:

    @staticmethod
    def _make_pdbt_line(
        count, atom_name, res_name, chain, res_num, coord, charge, atom_type, icode=""
    ):
        if atom_type is None:
            msg = "Can't write PDBT because atom_type is None.\n"
            msg += "Use PDBT atom typing (pdbt_atomtyper.assign_pdbt_types) before writing."
            raise ValueError(msg)
        x, y, z = float(coord[0]), float(coord[1]), float(coord[2])
        line = "ATOM  %5d %-4s %-3s %s%4d%s   %8.3f%8.3f%8.3f  0.00  0.00    %+5.3f %.2s" % (
            count,
            atom_name,
            res_name,
            chain if chain else " ",
            res_num,
            icode if icode else " ",
            x,
            y,
            z,
            charge,
            atom_type,
        )
        return line

    @classmethod
    def _make_pdbt_line_from_molsetup(cls, setup, atom_idx, count):
        pdbinfo = setup.get_pdbinfo(atom_idx)
        if pdbinfo is None:
            pdbinfo = pdbutils.PDBAtomInfo("", "", 0, "")
        atom_name, res_name, res_num, chain = cls._get_pdbinfo_fitting_pdb_chars(pdbinfo)
        if not atom_name.strip():
            atom_name = cls._generate_atom_name(setup, atom_idx)
        coord = setup.get_coord(atom_idx)
        atom_type = setup.get_atom_type(atom_idx)
        charge = setup.get_charge(atom_idx)
        return cls._make_pdbt_line(
            count, atom_name, res_name, chain, res_num, coord, charge, atom_type
        )

    @staticmethod
    def _generate_atom_name(setup, atom_idx):
        atomic_num = setup.get_atomic_num(atom_idx)
        symbol = _ELEMENT_SYMBOLS.get(atomic_num, "")
        if not symbol:
            try:
                atom = Chem.Atom(atomic_num)
                symbol = atom.GetSymbol()
            except Exception:
                symbol = "X"
        return " %-3s" % symbol

    @staticmethod
    def _get_pdbinfo_fitting_pdb_chars(pdbinfo):
        atom_name = pdbinfo.name
        res_name = pdbinfo.resName
        res_num = pdbinfo.resNum
        chain = pdbinfo.chain
        if len(atom_name) > 4:
            atom_name = atom_name[0:4]
        if len(res_name) > 3:
            res_name = res_name[0:3]
        if res_num > 9999:
            res_num = res_num % 10000
        if len(chain) > 1:
            chain = chain[0:1]
        return atom_name, res_name, res_num, chain

    @staticmethod
    def _is_molsetup_ok(setup, bad_charge_ok):
        success = True
        error_msg = ""
        if len(setup.restraints):
            error_msg = "molsetup has restraints but these can't be written to PDBT"
            success = False
        for atom in setup.atoms:
            if atom.is_ignore:
                continue
            if atom.atom_type is None:
                error_msg += "atom number %d has None type, mol name: %s\n" % (
                    atom.index,
                    setup.get_mol_name(),
                )
                success = False
        for atom in setup.atoms:
            if atom.is_ignore:
                continue
            if atom.atom_type is None:
                error_msg += "atom number %d has None type, mol name: %s\n" % (
                    atom.index,
                    setup.get_mol_name(),
                )
                success = False
            c = atom.charge
            if not bad_charge_ok and (
                type(c) not in (float, int) or math.isnan(c) or math.isinf(c)
            ):
                error_msg += (
                    "atom number %d has non finite charge, mol name: %s, charge: %s, type: %s\n"
                    % (atom.index, setup.get_mol_name(), str(c), type(c))
                )
                success = False
        return success, error_msg

    @classmethod
    def _collect_tree(cls, setup):
        """Walk flexibility model and collect atoms in tree order with structure entries."""
        data = {
            "visited": [],
            "atoms": [],
            "entries": [],
        }
        cls._collect_recursive(setup, setup.flexibility_model["root"], data, first=True)
        return data

    @classmethod
    def _collect_recursive(cls, setup, node, data, edge_start=0, first=False):
        if first:
            data["entries"].append({"type": "ROOT"})
            members = sorted(setup.flexibility_model["rigid_body_members"][node])
        else:
            members = setup.flexibility_model["rigid_body_members"][node][:]
            members.remove(edge_start)
            members = [edge_start] + members

        for member in members:
            if setup.get_is_ignore(member):
                continue
            if setup.atoms[member].is_pseudo_atom:
                continue
            data["atoms"].append(member)
            data["entries"].append({"type": "ATOM", "pos": len(data["atoms"]) - 1})

        if first:
            data["entries"].append({"type": "ENDROOT"})

        data["visited"].append(node)

        for neigh in setup.flexibility_model["rigid_body_graph"][node]:
            if neigh in data["visited"]:
                continue
            begin_idx, next_idx = setup.flexibility_model["rigid_body_connectivity"][node, neigh]
            if setup.get_is_ignore(begin_idx) or setup.get_is_ignore(next_idx):
                continue
            if setup.atoms[begin_idx].is_pseudo_atom or setup.atoms[next_idx].is_pseudo_atom:
                continue

            branch_start_pos = len(data["atoms"])
            data["entries"].append({
                "type": "BRANCH",
                "begin_idx": begin_idx,
                "branch_start_pos": branch_start_pos,
            })
            cls._collect_recursive(setup, neigh, data, edge_start=next_idx)
            data["entries"].append({
                "type": "ENDBRANCH",
                "begin_idx": begin_idx,
                "branch_start_pos": branch_start_pos,
            })

    @classmethod
    def write_string(cls, setup, bad_charge_ok=False):
        success, error_msg = cls._is_molsetup_ok(setup, bad_charge_ok)
        if not success:
            return "", success, error_msg

        tree = cls._collect_tree(setup)

        # Filter: remove non-HD H atoms
        keep_positions = set()
        pos_to_new_serial = {}
        serial = 1
        for pos, atom_idx in enumerate(tree["atoms"]):
            is_filtered = setup.get_atomic_num(atom_idx) == 1 and setup.get_atom_type(atom_idx) != "HD"
            if not is_filtered:
                keep_positions.add(pos)
                pos_to_new_serial[pos] = serial
                serial += 1

        # Calculate torsion info
        active_torsions = 0
        torsion_info = []
        for entry in tree["entries"]:
            if entry["type"] == "BRANCH":
                begin_idx = entry["begin_idx"]
                branch_start_pos = entry["branch_start_pos"]
                begin_pos = tree["atoms"].index(begin_idx)
                if begin_pos not in keep_positions:
                    continue
                first_kept = None
                for pos in sorted(keep_positions):
                    if pos >= branch_start_pos:
                        first_kept = pos
                        break
                if first_kept is None:
                    continue
                begin_serial = pos_to_new_serial[begin_pos]
                first_serial = pos_to_new_serial[first_kept]
                active_torsions += 1
                begin_elem = _ELEMENT_SYMBOLS.get(setup.get_atomic_num(begin_idx), "X")
                next_elem = _ELEMENT_SYMBOLS.get(setup.get_atomic_num(tree["atoms"][first_kept]), "X")
                torsion_info.append((active_torsions, begin_serial, begin_elem, first_serial, next_elem))

        buffer = []
        buffer.append("REMARK  Name = %s" % setup.get_mol_name())
        if active_torsions > 0:
            buffer.append("REMARK  %d active torsions:" % active_torsions)
            buffer.append("REMARK  status: ('A' for Active; 'I' for Inactive)")
            for idx, begin_serial, begin_elem, first_serial, next_elem in torsion_info:
                buffer.append("REMARK  %3d  A    between atoms: %s_%d  and  %s_%d" % (
                    idx, begin_elem, begin_serial, next_elem, first_serial
                ))
        buffer.append("REMARK                            x       y       z     vdW  Elec       q    Type")
        buffer.append("REMARK                         _______ _______ _______ _____ _____    ______ ____")

        for entry in tree["entries"]:
            if entry["type"] == "ROOT":
                buffer.append("ROOT")
            elif entry["type"] == "ENDROOT":
                buffer.append("ENDROOT")
            elif entry["type"] == "ATOM":
                pos = entry["pos"]
                if pos in pos_to_new_serial:
                    new_serial = pos_to_new_serial[pos]
                    atom_idx = tree["atoms"][pos]
                    line = cls._make_pdbt_line_from_molsetup(setup, atom_idx, new_serial)
                    buffer.append(line)
            elif entry["type"] == "BRANCH":
                begin_idx = entry["begin_idx"]
                branch_start_pos = entry["branch_start_pos"]
                begin_pos = tree["atoms"].index(begin_idx)
                if begin_pos not in keep_positions:
                    continue
                first_kept = None
                for pos in sorted(keep_positions):
                    if pos >= branch_start_pos:
                        first_kept = pos
                        break
                if first_kept is None:
                    continue
                begin_serial = pos_to_new_serial[begin_pos]
                first_serial = pos_to_new_serial[first_kept]
                buffer.append("BRANCH %3d %3d" % (begin_serial, first_serial))
            elif entry["type"] == "ENDBRANCH":
                begin_idx = entry["begin_idx"]
                branch_start_pos = entry["branch_start_pos"]
                begin_pos = tree["atoms"].index(begin_idx)
                if begin_pos not in keep_positions:
                    continue
                first_kept = None
                for pos in sorted(keep_positions):
                    if pos >= branch_start_pos:
                        first_kept = pos
                        break
                if first_kept is None:
                    continue
                begin_serial = pos_to_new_serial[begin_pos]
                first_serial = pos_to_new_serial[first_kept]
                buffer.append("ENDBRANCH %3d %3d" % (begin_serial, first_serial))

        buffer.append("TORSDOF %d" % active_torsions)

        pdbt_string = "\n".join(buffer) + "\n"
        return pdbt_string, success, error_msg
