# 🧬 Installing and Using Open Babel

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


# ⌨️ smiles2xyz.py - SMILES to 3D Structure Converter

This script serves as the starting point for the STEPs automated parameterization pipeline. It reads SMILES strings from a CSV file and converts them into 3D geometries (`.xyz`) for quantum mechanics calculations, while simultaneously generating 2D images for visual inspection.

## Features

1.  **3D Coordinate Generation (`.xyz`):** Uses OpenBabel to convert 1D SMILES strings into 3D structures required for **xTB/ORCA** optimization [Source 4, 5].
2.  **Visual Inspection (`.png`):** Uses RDKit to generate 2D structural diagrams labeled with atom indices [Source 6].
3.  **Connectivity Data (`.sdf`):** Preserves bond information and atom ordering [Source 4].
4.  **Atom Mapping (`.csv`):** Creates a mapping file linking atom indices to element types [Source 8].

## Prerequisites

This script requires **Python 3**, **OpenBabel**, and **RDKit**. We recommend installing them via Conda:

```bash
conda install -c conda-forge openbabel rdkit pandas
```
*Note: The script explicitly checks for the `obabel` executable in your system PATH [Source 2].*

## Quick Start

For the STEPs workflow (compatible with **Amber/mdgx**), use the following command to generate 3D coordinates and use **1-based indexing**:

```bash
python smiles2xyz.py --csv data.csv --gen3d --index-base 1
```

### Why these flags?
*   `--gen3d`: **Crucial.** Tells OpenBabel to generate actual 3D coordinates. Without this, your molecule will be flat, which is bad for xTB optimization [Source 4].
*   `--index-base 1`: **Crucial.** Sets atom numbering to start at 1 (instead of the default 0). This ensures the atom IDs in the PNG images match the IDs used in Amber `tleap`, `mdgx`, and Gaussian inputs [Source 7, 9].

## Usage & Arguments

```bash
python smiles2xyz.py [arguments]
```

| Argument | Default | Description |
| :--- | :--- | :--- |
| `--csv` | *(Required)* | Path to the input CSV file containing SMILES. |
| `--smiles-col` | `Capped_SMILES` | The column name in the CSV holding the SMILES string [Source 11]. |
| `--gen3d` | `False` | **Recommended.** If set, generates 3D coordinates using OpenBabel. |
| `--index-base` | `0` | Starting number for atom indices (`0` or `1`). **Set to `1` for Amber/mdgx.** |
| `--outdir` | `out` | Directory where output files will be saved. |
| `--prefix` | `STEPs` | Prefix for output filenames (e.g., `STEPs_0001.xyz`). |
| `--png-size` | `450` | Size of the output PNG images in pixels. |

## Output Structure

The script creates an output directory (default: `out/`) with the following subfolders [Source 16]:

*   **`out/xyz/`**: Contains `.xyz` files.
    *   *Use:* Input for **xTB** or **ORCA** geometry optimization.
*   **`out/png/`**: Contains `.png` images.
    *   *Use:* Visual check of the structure and atom indices (for setting up `mdgx` constraints).
*   **`out/sdf/`**: Contains `.sdf` files.
    *   *Use:* Intermediate file preserving connectivity.
*   **`out/atommap/`**: Contains `_atommap.csv` files.
    *   *Use:* Detailed list of atoms and their corresponding indices.

## FAQ

**Q: Why do I need to set `--index-base 1`?**
A: Python and RDKit use 0-based indexing by default [Source 11]. However, molecular dynamics software like Amber (`mdgx`) and Gaussian use 1-based indexing. Setting this flag ensures the numbers you see on the generated PNG images exactly match the numbers you need to write in your MD input files.

**Q: My input CSV column isn't named "Capped_SMILES".**
A: You can specify your column name using the argument: `--smiles-col "YourColumnName"`.

**Q: I get a "Cannot find 'obabel'" error.**
A: The script cannot find the OpenBabel executable [Source 2]. Ensure you have installed OpenBabel and activated your conda environment.
