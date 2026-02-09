# mdgx Angle/Dihedral 2D Map Generator (RDKit + PIL)

This script reads an **mdgx input file** (e.g., `genconf.in` / `generate.in`), extracts each `GridSample` and `RandomSample` line, and generates **one 2D structure PNG per config line**.

Each PNG shows:

* The **same reference molecule** drawn as a 2D structure
* Only the **atoms involved in that config line** labeled by **atom name**
* Bonds along the specified atom sequence highlighted
* A legend with sampling type, atom sequence, range, and `Krst`

It is designed for quickly verifying:

* Which atom names are used for each **dihedral** (`GridSample`, 4 atoms)
* Which atom names are used for each **angle** (`RandomSample`, 3 atoms)

---

## Requirements

### Python packages

* `rdkit`
* `pillow` (PIL)

Recommended installation via conda:

```bash
conda create -n mdgx_draw python=3.10 -y
conda activate mdgx_draw
conda install -c conda-forge rdkit pillow -y
```

Quick check:

```bash
python -c "from rdkit import Chem; from PIL import Image; print('OK')"
```

---

## Input files

You need **two** inputs:

1. **mdgx input file** (e.g., `genconf.in` or `generate.in`)

   * Must contain lines like:

     * `GridSample :1@C2 :1@N1 :1@C10 :1@C11 { -180.0 180.0 } Krst 64.0`
     * `RandomSample :1@C2 :1@N1 :1@C3 { 105.0 135.0 } Krst 256.0`

2. **Reference structure file** (`.mol2` or `.pdb`)

   * Used to draw the 2D structure and map atom names → atom indices.

> The script **does not need conformer files**.
> It only uses the reference structure as a drawing template.

---


## Usage

### Basic (auto-select reference in current directory)

```bash
python mdgx_atommap.py gen_coor.mdgx conf_png/
```

The script will pick a reference structure in the current directory with this priority:

1. first **non**-`Conf*.mol2`
2. first **non**-`Conf*.pdb`

### Explicit reference (recommended)

```bash
python mdgx_atommap.py gen_coor.mdgx conf_png/ R20.mol2
```

### Output

If your mdgx config contains 7 sampling lines, you will get 7 PNGs:

```text
angle_maps/
  01_GridSample_C2-N1-C10-C11.png
  02_GridSample_N1-C10-C11-N2.png
  ...
  06_RandomSample_C2-N1-C3.png
  07_RandomSample_C10-N1-C3.png
```

---

## Running mdgx to generate conformers

After preparing `generate.in`, conformer generation is performed using **mdgx**, which is part of **AmberTools**.

### Prerequisites

* AmberTools installed and available in your environment
  (verify with: `which mdgx`)
* A topology and coordinate file for the reference molecule:

  * `R20.top`
  * `R20.crd`

Typical directory layout:

```text
R20_example/
  generate.in
  R20.top
  R20.crd
  R20.mol2              # reference structure (for visualization only)
  mdgx_angles_2Dmaps.py
```

---

### Minimal `generate.in` structure

A typical `generate.in` used for conformer generation looks like:

```text
&files
  -p    ./R20.top
  -c    ./R20.crd
  -o    GenConformers.out
&end

&configs
  GridSample    :1@C2 :1@N1 :1@C10 :1@C11 { -180.0 180.0 } Krst 64.0
  RandomSample  :1@C2 :1@N1 :1@C3        { 105.0 135.0 }  Krst 256.0

  combine 1 2
  count 3000
  verbose 1

  qmlev    'B3LYP',
  basis    '6-31G**',

  outbase  'Conf', 'Conf'
  write    'pdb',  'orca'
  outsuff  'pdb',  'orca'
&end
```

Only the `&configs` block controls **how conformers are sampled**.
The Python visualization script reads **only this block**.

---

### Running mdgx

From the directory containing `generate.in`:

```bash
mdgx -i gen_coor.mdgx
```

mdgx will:

1. Read the reference topology (`.top`) and coordinates (`.crd`)
2. Apply the sampling definitions in `&configs`
3. Generate a set of restrained conformers
4. Write conformers to individual files:

   ```text
   Conf000001.pdb
   Conf000002.pdb
   ...
   ```

The run log is written to:

```text
GenConformers.out
```

---

### Notes on output files

* Each `Conf*.pdb` corresponds to **one conformer**
* The numbering reflects the order in which conformers are generated
* The PDB files may later be used for:

  * QM single-point calculations
  * RESP charge derivation
  * Force-field parameter fitting

⚠️ **Important**
The Python script **does not read these conformer files**.
It uses only the *reference structure* and *generate.in* to visualize the sampling design.

---




### Summary of roles in the workflow

```text
R20.mol2
   ↓
gen_coor.mdgx  (defines WHAT to sample)
   ↓
mdgx         (generates conformers)
   ↓
Conf*.pdb    (conformer ensemble)

mdgx_atommap.py
   ↳ visualizes sampling definitions only
```

---

## Is this code only for MOL2?

**No.** It supports both:

* ✅ **MOL2** (`.mol2`) — recommended

  * Atom names are read from `_TriposAtomName`
  * Usually matches `:1@C10` style naming reliably

* ✅ **PDB** (`.pdb` / `.ent`)

  * Atom names are read from PDB residue atom name field (`ATOM` name column)

However, **MOL2 is strongly recommended** because:

* Some PDB-writing tools may rename atoms or change formatting
* MOL2 atom naming tends to stay consistent with force-field workflows

### Important limitation (applies to both)

If your reference file’s atom names do not match the mdgx input (e.g., reference uses `C1` but mdgx uses `C10`), the script will warn:

```text
[WARN] ... missing atom names: ['C10', ...]
```

In that case, fix the reference naming (or use the correct reference file).

---

## Customization

Inside `draw_one()` you can adjust:

* Canvas size:

  ```python
  size = (2200, 1600)
  ```

* 2D layout scale (bigger = more spread out):

  ```python
  AllChem.Compute2DCoords(mol, bondLength=2.2)
  ```

* Label placement (PIL overlay):

  ```python
  off = 75
  extra_y = [-40, 40, -40, 40]
  ```

If labels go off-screen, reduce `off` or use a larger canvas.

---

## Troubleshooting

### 1) “No GridSample/RandomSample lines found”

Make sure your mdgx input contains literal lines starting with:

* `GridSample`
* `RandomSample`

### 2) Missing atom names warnings

Your mdgx input atom names don’t exist in the reference structure naming.

* Use a different reference file (prefer `.mol2`)
* Or rename atoms in the reference to match mdgx input

### 3) RDKit not found

Install with conda-forge:

```bash
conda install -c conda-forge rdkit pillow -y
```

---
