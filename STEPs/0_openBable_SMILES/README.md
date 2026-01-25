## Installing and Using Open Babel

**Open Babel** is used for molecular file format conversion and initial structure generation (e.g., SMILES → XYZ / PDB), which serves as the starting point for subsequent xTB calculations.

---

## 1. Install Open Babel via Conda

Open Babel can be installed directly from `conda-forge` and does not require root privileges.

Activate the desired Conda environment (e.g., the same environment used for xTB):

```bash
conda activate xtb
```

Install Open Babel:

```bash
conda install openbabel
```

Verify installation:

```bash
obabel -V
```

A successful installation will print the Open Babel version.

---

## 2. Generate Initial Structures from SMILES

### 2.1 SMILES → XYZ (with 3D coordinates)

```bash
obabel -:"CCN(CC)CC" -O molecule.xyz --gen3d
```

* `-:"SMILES"` : inline SMILES input
* `-O molecule.xyz` : output file
* `--gen3d` : generate 3D coordinates

---

### 2.2 SMILES → PDB

```bash
obabel -:"CCN(CC)CC" -O molecule.pdb --gen3d
```

This is useful when downstream tools (e.g., AmberTools, VMD) require PDB input.

---

## 3. File Format Conversion

Open Babel supports a wide range of molecular formats.

### 3.1 XYZ → PDB

```bash
obabel molecule.xyz -O molecule.pdb
```

### 3.2 MOL2 → XYZ

```bash
obabel molecule.mol2 -O molecule.xyz
```

---

## 4. Add or Modify Hydrogen Atoms

### Add hydrogens (based on valence rules)

```bash
obabel molecule.xyz -O molecule_H.xyz -h
```

### Remove hydrogens

```bash
obabel molecule.xyz -O molecule_noH.xyz -d
```

---

## 5. Geometry Preprocessing (Optional)

Open Babel can perform a quick force-field-based cleanup before quantum calculations.

```bash
obabel molecule.xyz -O molecule_preopt.xyz --minimize --ff MMFF94
```

*Note:* This step is optional. In this workflow, final geometries are optimized using xTB.

---

## 6. Typical Workflow Integration

A common preprocessing pipeline is:

```text
SMILES
  ↓ (Open Babel: --gen3d)
XYZ / PDB
  ↓ (xTB: GFN2-xTB)
Optimized structure
```

Example:

```bash
obabel -:"CCN(CC)CC" -O molecule.xyz --gen3d
xtb molecule.xyz --gfn 2
```

---

## Notes

* Open Babel is used **only for structure generation and format conversion**.
* Energetic evaluations and geometry optimizations are performed using *GFN2-xTB*.
* The initial geometry from Open Babel does not affect final results after xTB optimization.

---

## Minimal Version (for quick README sections)

```markdown
Open Babel is used to generate initial 3D structures from SMILES:

obabel -:"SMILES" -O molecule.xyz --gen3d
```

