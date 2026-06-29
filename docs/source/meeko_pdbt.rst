
Meeko-PDBT: AutoDock PDBT Output
==================================

Overview
--------

:command:`meeko-pdbt` is a modified version of Meeko that produces
AutoDock PDBT (Protein Data Bank, Partial Charge (Q), & Atom Type (T)) output
files, mimicking the behavior of OpenBabel 25.07.

The PDBT format is used by AutoDock-family docking tools and includes:

- REMARK header with molecule name, active torsion listing, and column headers
- ROOT/BRANCH/ENDBRANCH/ENDROOT torsion tree structure
- ATOM lines with coordinates, partial charges, and two-letter Vinardo2 atom types
- TORSDOF footer with total torsional degrees of freedom

Quickstart
----------

Command line: prepare a PDBT ligand from an SDF
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``meeko-pdbt`` command is installed automatically with Meeko. It accepts
the same arguments as ``mk_prepare_ligand.py`` (it is a thin wrapper that
adds the PDBT writer), so any option that works for the PDBQT ligand
preparation also works here.

.. code-block:: bash

   # Print the help
   meeko-pdbt -h

   # Write a PDBT file (extension defaults to .pdbt)
   meeko-pdbt -i ligand.sdf -o ligand.pdbt

   # Or to stdout (useful for piping into docking tools)
   meeko-pdbt -i ligand.sdf -o - > ligand.pdbt

The input SDF must contain 3D coordinates and explicit hydrogens (or
none, in which case RDKit's ``AddHs`` is called). For best results,
match obabel's ``-p7`` protonation model by pre-titrating the ligand
with your cheminformatics toolkit of choice.

Python API: write a PDBT string
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For finer control (e.g. to embed the PDBT string directly into another
Python pipeline, or to add custom atom typing), use the Python API.
``PDBTWriterLegacy`` is exported from the top-level ``meeko`` package:

.. code-block:: python

   from rdkit import Chem
   from rdkit.Chem import AllChem
   from meeko import MoleculePreparation
   from meeko import PDBTWriterLegacy
   from meeko.pdbt_atomtyper import assign_pdbt_types

   # Build and embed a 3D molecule with explicit Hs
   mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")  # aspirin
   mol = Chem.AddHs(mol)
   AllChem.EmbedMolecule(mol)
   AllChem.MMFFOptimizeMolecule(mol)

   # Prepare (same as for PDBQT)
   preparator = MoleculePreparation()
   molsetups = preparator.prepare(mol)
   molsetup = molsetups[0]

   # Assign Vinardo2 / PDBT atom types
   assign_pdbt_types(molsetup, mol)

   # Write the PDBT string
   pdbt_string, success, error_msg = PDBTWriterLegacy.write_string(molsetup)
   if not success:
       raise RuntimeError(error_msg)

   with open("aspirin.pdbt", "w") as f:
       f.write(pdbt_string)

The ``assign_pdbt_types`` step is essential: without it the molsetup
still holds AD4 types (e.g. ``A``, ``C``, ``OA``) and the writer will
refuse to emit a PDBT line.

Command-line help
^^^^^^^^^^^^^^^^^

.. code-block:: text

   usage: mk_prepare_pdbt_ligand.py [-h] [-v] -i INPUT_MOLECULE_FILENAME
                                    [-o OUTPUT_PDBT_FILENAME] [-]
                                    [-c CONFIG_FILE] [--rigid_macrocycles]
                                    [--keep_chorded_rings]
                                    [--keep_equivalent_rings]
                                    [--min_ring_size MIN_RING_SIZE]
                                    [-r SMARTS] [-b i j [i j ...]] [-a]
                                    [--double_bond_penalty DOUBLE_BOND_PENALTY]
                                    [--charge_model {gasteiger,zero,read}]
                                    [--charge_atom_prop CHARGE_ATOM_PROP]
                                    [--bad_charge_ok] [--rename_atoms]

   options:
     -h, --help            show this help message and exit
     -v, --verbose         print information about molecule setup

   Input/Output:
     -i INPUT_MOLECULE_FILENAME, --mol INPUT_MOLECULE_FILENAME
                           molecule file (MOL2, SDF, ...)
     -o OUTPUT_PDBT_FILENAME, --out OUTPUT_PDBT_FILENAME
                           output pdbt filename. Single molecule input only.
     -, --                 do not write file, redirect output to STDOUT.

   Molecule preparation:
     -c CONFIG_FILE, --config_file CONFIG_FILE
                           configure MoleculePreparation from JSON file.
     --rigid_macrocycles   keep macrocycles rigid in input conformation
     --keep_chorded_rings  return all rings from exhaustive perception
     --keep_equivalent_rings
                           equivalent rings have the same size and neighbors
     --min_ring_size MIN_RING_SIZE
                           min nr of atoms in ring for opening
     -r SMARTS, --rigidify_bonds_smarts SMARTS
                           SMARTS patterns to rigidify bonds
     -b i j [i j ...], --rigidify_bonds_indices i j [i j ...]
                           indices of two atoms (in the SMARTS) that define a
                           bond (start at 1)
     -a, --flexible_amides
                           allow amide bonds to rotate and be non-planar, which
                           is bad
     --double_bond_penalty DOUBLE_BOND_PENALTY
                           penalty > 100 prevents breaking double bonds
     --charge_model {gasteiger,zero,read}
                           default is 'gasteiger', 'zero' sets all zeros
     --charge_atom_prop CHARGE_ATOM_PROP
                           set atom partial charges from an RDKit atom property
                           based on the input file.
     --bad_charge_ok       NaN and Inf charges allowed in PDBT
     --rename_atoms        rename atoms: new name is original name + (1-based)
                           index

Receptor PDBT output
^^^^^^^^^^^^^^^^^^^^

The :command:`meeko-pdbt-rec` CLI writes a receptor PDBT from a PDB or
mmCIF file. Like the ligand CLI, it is intentionally simple: no
flexres, no reactive residues, no GPF/box output. It exists
independently from :command:`mk_prepare_receptor.py` (which only
emits PDBQT).

Command line: prepare a PDBT receptor from a PDB file
+++++++++++++++++++++++++++++++++++++++++++++++++++++++

.. code-block:: bash

   # Write a receptor PDBT (extension defaults to .pdbt)
   meeko-pdbt-rec -i receptor.pdb -o receptor.pdbt

   # Read mmCIF (requires ProDy) and write to stdout
   meeko-pdbt-rec -i receptor.cif -o -

   # Handle alternate locations, delete bad residues, etc.
   meeko-pdbt-rec -i receptor.pdb -o receptor.pdbt \
       --default_altloc A -x

Command-line help
+++++++++++++++++

.. code-block:: text

   usage: mk_prepare_pdbt_receptor.py [-h] -i INPUT_FILENAME
                                       [-o OUTPUT_FILENAME] [-]
                                       [--default_altloc DEFAULT_ALTLOC]
                                       [--wanted_altloc WANTED_ALTLOC]
                                       [-n SET_TEMPLATE]
                                       [-d DELETE_RESIDUES]
                                       [-b BLUNT_ENDS]
                                       [--config_file CONFIG_FILE]
                                       [--add_templates ADD_TEMPLATES]
                                       [-x] [--forgive_extra_bonds]
                                       [--bad_charge_ok]

   Input/Output:
     -i INPUT_FILENAME, --input INPUT_FILENAME
                           receptor file (PDB or mmCIF; mmCIF requires ProDy)
     -o OUTPUT_FILENAME, --output OUTPUT_FILENAME
                           output PDBT filename. Default: <input>.pdbt.
                           Use '-' for stdout.
     -, --stdout           do not write file, redirect output to STDOUT

   Receptor perception:
     --default_altloc DEFAULT_ALTLOC
                           default alternate location (e.g. 'A')
     --wanted_altloc WANTED_ALTLOC
                           require altloc for specific residues,
                           e.g. :5=B,B:17=A
     -n SET_TEMPLATE, --set_template SET_TEMPLATE
                           override residue template, e.g. A:5,7=CYX
     -d DELETE_RESIDUES, --delete_residues DELETE_RESIDUES
                           delete residues by chain:number,
                           e.g. A:350,B:15,16,17
     -b BLUNT_ENDS, --blunt_ends BLUNT_ENDS
                           mark chain ends as blunt, e.g. A:123=2,A:1=0
     --config_file CONFIG_FILE
                           JSON file with extra MoleculePreparation
                           settings
     --add_templates ADD_TEMPLATES
                           additional residue template (JSON file or
                           'resname:file.sdf'); repeatable
     -x, --delete_bad_res   delete residues that don't match templates
                            instead of raising
     --forgive_extra_bonds  allow processing structures with excess
                            bonds (use with care)
     --bad_charge_ok        allow NaN/Inf charges in the output

Python API: write a receptor PDBT
+++++++++++++++++++++++++++++++++

.. code-block:: python

   from meeko import MoleculePreparation, Polymer, ResidueChemTemplates
   from meeko import PDBTWriterLegacy
   from meeko.pdbt_atomtyper import assign_pdbt_types_from_pdbinfo

   # Read a receptor PDB
   with open("receptor.pdb") as f:
       pdb_string = f.read()

   # Build the polymer
   templates = ResidueChemTemplates.create_from_defaults()
   mk_prep = MoleculePreparation()
   polymer = Polymer.from_pdb_string(pdb_string, templates, mk_prep)

   # Assign PDBT types using the hardcoded protein lookup
   # (with chemical-rules fallback for non-standard residues)
   assign_pdbt_types_from_pdbinfo(polymer)

   # Write a single flat receptor PDBT
   pdbt_string, success, error_msg = PDBTWriterLegacy.write_from_polymer(polymer)
   if not success:
       raise RuntimeError(error_msg)

   with open("receptor.pdbt", "w") as f:
       f.write(pdbt_string)

The output mirrors what OpenBabel-25-07 produces for receptors: one
``ROOT/ENDROOT`` containing all heavy atoms + HD hydrogens, no torsion
tree (``TORSDOF 0``). Atoms marked ``is_ignore`` (e.g. non-polar Hs)
and macrocycle glue atoms (``is_pseudo_atom``) are skipped.

Implementation
--------------

Files Created/Modified
^^^^^^^^^^^^^^^^^^^^^^

+--------------------------------------+----------------------------------------------------+
| File                                 | Description                                        |
+======================================+====================================================+
| ``meeko/pdbt_atomtyper.py``          | Python implementation of PDBT/Vinardo2 atom typing |
+--------------------------------------+----------------------------------------------------+
| ``meeko/pdbt_writer.py``            | PDBT output writer with REMARK header and torsion  |
|                                      | tree                                               |
+--------------------------------------+----------------------------------------------------+
| ``meeko/cli/mk_prepare_pdbt_ligand.py`` | CLI entry point for ``meeko-pdbt``             |
+--------------------------------------+----------------------------------------------------+
| ``meeko/__init__.py``               | Exports ``PDBTWriterLegacy``                       |
+--------------------------------------+----------------------------------------------------+
| ``setup.py``                        | Registers ``meeko-pdbt`` console_scripts entry     |
|                                      | point                                              |
+--------------------------------------+----------------------------------------------------+

Atom Typing
^^^^^^^^^^^

The PDBT atom typing (:file:`pdbt_atomtyper.py`) implements the Vinardo2
two-letter atom type scheme from OpenBabel 25.07's :file:`pdbtformat.cpp`.
It uses RDKit to determine atom properties and maps them to the same type
classes used by OpenBabel.

**Atom Type Classes:**

.. list-table::
   :header-rows: 1

   * - Type
     - Description
   * - A0
     - Aromatic carbon, 0 hetero neighbors
   * - A1
     - Aromatic carbon, 1 hetero neighbor
   * - C1
     - Aliphatic carbon, hvy_deg 1 (sp3 C with 1 heavy neighbor)
   * - C2
     - Aliphatic carbon, hvy_deg 2 (sp3 C with 2 heavy neighbors)
   * - C3
     - Aliphatic carbon, hvy_deg 3 (sp3 C with 3 heavy neighbors)
   * - Np
     - sp2 nitrogen, acceptor, in ring, no hetero neighbor
   * - Nf
     - sp3 nitrogen, donor
   * - Ng
     - sp2 nitrogen, guanidinium-like
   * - Nh
     - sp2 nitrogen, donor
   * - Ni
     - sp3 nitrogen, not donor, not acceptor
   * - Nj
     - sp3 nitrogen, not donor, acceptor
   * - Nk
     - sp3 nitrogen, donor, acceptor
   * - Nl
     - sp3 nitrogen, donor, not acceptor
   * - Nu
     - sp2 nitrogen, not donor, not acceptor
   * - Oa
     - sp3 oxygen, acceptor (hydroxyl)
   * - Ob
     - sp2 oxygen, donor+acceptor (carboxyl)
   * - Oc
     - sp2 oxygen, acceptor (carbonyl)
   * - Oe
     - sp2 oxygen, acceptor (ester carbonyl)
   * - Of
     - sp2 oxygen, other
   * - Og
     - sp2 oxygen, other
   * - Ok
     - sp2 oxygen, other (nitro)
   * - Ol
     - sp3 oxygen, other
   * - Oo
     - sp2 oxygen, other (ether-like)
   * - S3
     - sp3 sulfur
   * - S2
     - sp2 sulfur (thiophene)
   * - F
     - Fluorine
   * - Cl
     - Chlorine
   * - Br
     - Bromine
   * - I
     - Iodine
   * - HD
     - Polar hydrogen (on N, O, S)

Key Design Decisions
^^^^^^^^^^^^^^^^^^^^

1. **Pure Python implementation**: PDBT atom typing is implemented in pure
   Python with SMARTS-like logic derived from the C++ functions, rather than
   using SMARTS pattern JSON files, for maximum fidelity to the OpenBabel
   logic.

2. **Charges preserved**: Unlike OpenBabel's PDBT output (which hardcodes
   charges to 0.000), ``meeko-pdbt`` outputs actual Gasteiger charges computed
   by Meeko, providing more useful information for scoring.

3. **Hydrogen filtering**: Non-polar hydrogens (non-HD type, i.e., carbon-bound
   H's) are filtered out of the output, matching OpenBabel's behavior. Polar
   hydrogens (HD type on N, O, S) are preserved.

4. **Atom names**: When PDB info is not available, atom names are auto-generated
   from element symbols.

Testing
-------

Small Molecule Verification Against OpenBabel 25.07
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Tested on a diverse set of functional groups:

.. list-table::
   :header-rows: 1

   * - Molecule
     - Atom Types
     - Match
   * - Aspirin
     - A0, A1, C3, HD, Ob, Oc, Og, Oo
     - Yes
   * - Paracetamol
     - A0, A1, C3, HD, Nf, Oa, Of
     - Yes
   * - Benzamidine
     - A0, A1, HD, Ng
     - Yes
   * - Sulfanilamide
     - A0, HD, Nd, Nj, Ol, S3
     - Yes
   * - Nitrobenzene
     - A0, Nu, Ok
     - Yes
   * - Thiophene
     - A0, S2
     - Yes
   * - Caffeine
     - A0, C3, Np, Oe
     - No\ [#caffeine]_

.. [#caffeine] Caffeine: ``Nu`` (obabel) vs ``Np`` (meeko). Caused by
   different aromaticity perception between RDKit and OpenBabel for one
   nitrogen in the purine system. Expected minor difference.

Large-Scale Validation: runs-n-poses Dataset
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The `runs-n-poses <https://github.com/forlilab/runs-n-poses>`_ dataset
provides ~1426 protein-ligand complexes for benchmarking. The
``vinardo_inputs/`` directory contains PDBT files prepared with OpenBabel
25.07 for both ligands and receptors. These serve as a reference for
validating ``meeko-pdbt`` at scale.

Ligand Comparison
~~~~~~~~~~~~~~~~~

For each system, the ground truth SDF (from
``ground_truth/<system_id>/ligand_files/``) is processed through both
pipelines::

   meeko-pdbt:     SDF → RDKit AddHs → meeko-pdbt → meeko.pdbt
   obabel-25-07:   SDF → obabel -i sdf -o pdbt -p7 → obabel.pdbt

Comparison checks:

1. **Atom type agreement**: per-atom PDBT type comparison (expect >95%
   agreement)
2. **Known differences**: aromaticity disagreements between RDKit and
   OpenBabel in edge cases (purines, sulfonamides, some heterocycles)
3. **Charge differences**: meeko-pdbt outputs actual Gasteiger charges;
   OpenBabel hardcodes 0.000 for all atoms
4. **TORSDOF counts**: may differ due to different flexibility model
   heuristics

Receptor (Protein) Comparison
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For each system, the ground truth receptor structure (from
``ground_truth/<system_id>/receptor.cif``) is processed through both
pipelines::

   meeko-pdbt:      CIF → mk_prepare_receptor (--default_altloc A) → PDBT
   obabel-25-07:    CIF → obabel -i cif -o pdbt -p7 -xr → obabel.pdbt

The ``-p7`` flag adds hydrogens at pH 7, and ``-xr`` outputs a rigid
molecule (no torsion tree, which is appropriate for receptors).

**Meeko-specific handling:**

.. list-table::
   :header-rows: 1

   * - Feature
     - meeko
     - obabel
   * - **AltLoc A**
     - Parses altloc column (col 17); defaults to ``A`` via
       ``--default_altloc A``; reports residues needing altloc selection
     - Ignores altloc; takes first altloc encountered
   * - **Charges**
     - Outputs actual Gasteiger charges
     - Hardcodes 0.000
   * - **Hydrogen addition**
     - Uses RDKit's ``AddHs()`` at neutral pH
     - Uses OpenBabel's ``-p7`` pH model
   * - **Atom typing**
     - Pure Python PDBT atom typer
     - Native C++ from ``pdbtformat.cpp``

Comparison checks:

1. **Per-residue atom type agreement** for mainchain and sidechain atoms
2. **Altloc resolution**: meeko's ``--default_altloc A`` vs obabel's implicit
   first-altloc
3. **pH 7 protonation**: subtle differences between RDKit and OpenBabel pH
   models (e.g., histidine tautomers, carboxylate vs carboxylic acid)
4. **Atom types**: same aromaticity differences as ligands, plus
   residue-specific typing (e.g., ``Nf`` for backbone amide nitrogens, ``Of``
   for backbone carbonyls)

Measured Agreement
~~~~~~~~~~~~~~~~~~

.. list-table:: 1426-ligand validation results
   :header-rows: 1

   * - Metric
     - Value
   * - Overall atom type match rate
      - **99.1%**
   * - Total atoms compared
      - 40,867
   * - Type mismatches (expected, aromaticity)
      - 361
   * - Ligands with meeko preparation errors
      - 2 / 1426

The big jumps came from imitating OpenBabel-25-07's H-bond detection:

- **SP2 N with 3 valences** (was 1044, now 0): OpenBabel-25-07's
   ``IsHbondAcceptor()`` returns false for an sp2 N with 3 valences
   (e.g. pyrrole N1, imidazole N1, amide N) because the lone pair is
   in the pi system. Meeko now applies the same guard. The
   ``Np→Nu``, ``Nr→Nu``, ``Nq→Nu`` categories all dropped to zero.

Remaining mismatches (~361) are:

- **pH 7 protonation differences** (170): tertiary amines that obabel
   protonates to NH+ at pH 7 but RDKit keeps as neutral N. Affects
   ``Ns→Ni`` (and a few ``Ni→Nf``).
- **Ester/amide detection differences** (~115): ``Og→Oh``, ``Of→Oh``,
   ``Oc→Of`` — the C-O-C(=O) ester and C-N-C(=O) amide patterns are
   detected slightly differently. Meeko's detection is stricter.
- **SDF→SER misinterpretation** (~30): for some ligands, obabel-25-07's
   PDBT output misclassifies the first N of the SDF as a SER backbone
   N, then re-labels nearby O atoms as SER backbone atoms and adds ~6
   synthetic SER atoms. meeko correctly treats them as ligand atoms.
- **Obabel quirks** (a handful): obabel is sometimes inconsistent
   with its own H-bond donor detection (P-O-H sometimes donor,
   sometimes not), gives ``S1`` to S atoms with no H neighbor, and
   calls α-amino N+ of amino acids "amide" (``Nf``) due to a
   permissive ``IsAmideNitrogen()`` check.
- **Hydrogen positioning**: RDKit and OpenBabel place polar hydrogens at
  slightly different positions (~0.05 Å), requiring separate coordinate
  matching tolerances (0.01 Å heavy atoms, 0.5 Å hydrogens)
- **Charges differ on all atoms** (by design): meeko outputs actual Gasteiger
  charges; obabel hardcodes 0.000 for all atoms

Running the Test
~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Point the runs-n-poses symlink to the ground truth dataset
   ln -sf /path/to/runs-n-poses-datasets/ground_truth runs-n-poses

   # Batch validation (1426 systems) using vinardo_inputs reference PDBTs:
   python test/test_pdbt_batch.py

   # Receptor comparison (20 systems with pre-computed references):
   python test/test_pdbt_receptor_compare.py

   # Ligand comparison (example for one system)
   obabel -i sdf runs-n-poses/5s9z__1__1.A_1.B__1.R/ligand_files/1.R.sdf \
          -o pdbt -p7 -O obabel_ligand.pdbt
   meeko-pdbt -i runs-n-poses/5s9z__1__1.A_1.B__1.R/ligand_files/1.R.sdf \
              -o meeko_ligand.pdbt

   # Receptor CLI usage:
   meeko-pdbt-rec -i receptor.pdb -o receptor.pdbt

   # Receptor comparison (example for one system, using CIF)
   obabel -i cif runs-n-poses/5s9z__1__1.A_1.B__1.R/receptor.cif \
          -o pdbt -p7 -xr -O obabel_receptor.pdbt

   # For systems with altloc variants, specify per-residue:
   mk_prepare_receptor --read_with_prody \
       runs-n-poses/5s9z__1__1.A_1.B__1.R/receptor.cif \
       --wanted_altloc ":42=B,A:17=A" --default_altloc A \
       --write_pdb meeko_receptor_prepped.pdb

Caveats
-------

1. **Input must have explicit hydrogens**: Meeko requires SDF input with
   explicit hydrogen atoms. Use RDKit's ``Chem.AddHs()`` or OpenBabel to add
   hydrogens before processing.

2. **Flexibility model differences**: The torsion tree structure (which branch
   connects where) may differ between ``meeko-pdbt`` and OpenBabel due to
   different root-selection algorithms. TORSDOF counts may also differ (meeko
   may count more or fewer rotatable bonds than obabel for the same molecule).

3. **Aromaticity**: RDKit and OpenBabel may perceive aromaticity differently
   in edge cases (e.g., purine ring systems), leading to different atom types
   for a small fraction of atoms.

Receptor Validation
~~~~~~~~~~~~~~~~~~~~

The ``meeko-pdbt-rec`` CLI is validated against the obabel-25-07 reference
receptors in ``runs-n-poses/vinardo_inputs/receptors/<system>/``. The
20 systems that ship with pre-computed references all match at 100%
or 99.8%:

.. list-table::
   :header-rows: 1

   * - System
     - Match rate
   * - 5s9y, 5sau, 5sav, 5saw, 5sb2, 5s9z
     - 100.0% (exact match)
   * - 5sdh
     - 100.0% (10843 atoms)
   * - 5sdu, 5sdy, 5se0, 5se2, 5se3, 5se5, 5se6, 5se8, 5se9, 5sea, 5seb, 5sec, 5see
     - 99.8% (4–5 mismatches per system)
   * - **Overall (20 systems, 65037 atoms)**
     - **99.9%**

Residue atom types come from a hardcoded lookup table (O(1) dict
access) matching the 20 standard amino acids and HIS variants
(HIS, HID, HIE, HIP) from obabel-25-07's ``pdbtformat.cpp``. For
non-standard residues or HETATMs, the chemical-rules path is used as
a fallback.

Run the validation with::

   python test/test_pdbt_receptor_compare.py

Or test a specific system::

   python test/test_pdbt_receptor_compare.py --system 5s9y__1__1.A__1.K
