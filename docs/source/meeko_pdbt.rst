
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

CLI Usage
^^^^^^^^^

.. code-block:: bash

   # Basic usage (output to stdout):
   meeko-pdbt -i input.sdf -

   # Output to file:
   meeko-pdbt -i input.sdf -o output.pdbt

   # With explicit hydrogens (required):
   meeko-pdbt -i input_with_Hs.sdf -o output.pdbt

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

Expected Agreement
~~~~~~~~~~~~~~~~~~

The comparison should show:

- **Atom types**: ~95%+ agreement for ligands, ~98%+ for receptors (receptors
  have fewer exotic heterocycles)
- **Aromaticity mismatches**: The primary source of atom type differences,
  occurring in ring systems where RDKit and OpenBabel disagree on aromaticity
  (e.g., purines, thiazoles, certain fused heterocycles)
- **Hydrogen counts**: Potential differences due to distinct pH 7 models for
  titratable groups
- **Altloc differences**: Minimal impact when meeko defaults to altloc A, but
  some systems may have residues with only altloc B (requiring manual
  ``--wanted_altloc``)

Running the Test
~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Point the runs-n-poses symlink to the ground truth dataset
   ln -sf /path/to/runs-n-poses-datasets/ground_truth runs-n-poses

   # Ligand comparison (example for one system)
   obabel -i sdf runs-n-poses/5s9z__1__1.A_1.B__1.R/ligand_files/1.R.sdf \
          -o pdbt -p7 -O obabel_ligand.pdbt
   meeko-pdbt -i runs-n-poses/5s9z__1__1.A_1.B__1.R/ligand_files/1.R.sdf \
              -o meeko_ligand.pdbt

   # Receptor comparison (example for one system, using CIF)
   obabel -i cif runs-n-poses/5s9z__1__1.A_1.B__1.R/receptor.cif \
          -o pdbt -p7 -xr -O obabel_receptor.pdbt
   mk_prepare_receptor --read_with_prody \
       runs-n-poses/5s9z__1__1.A_1.B__1.R/receptor.cif \
       --default_altloc A --write_pdb meeko_receptor_prepped.pdb

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
