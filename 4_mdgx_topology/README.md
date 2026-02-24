# calculating single point energy with every conformers in Gaussian

## 1. create dihedral conformers
### 1.1 print mol2 2D structure figure to assign the atoms for dihedral

```bash
python mol2_atomname.py -i /path/to/s01.mol2 -o s01_names.png --hide-h-labels
```
- create png file without hydrogen
```bash
python mol2_atomname.py -i s01.mol2 -o s01_names.svg
```
- or svg file
---

### 1.3 run mdgx to create different conformers
example for mdgx input file
```bash
&files
  -p    ../../3_antechamber/S01/s01.top
  -c    ../../3_antechamber/S01/s01.crd
  -o    GenConformers.out
&end

&configs
  GridSample :1@C1 :1@N :1@C8 :1@C9 { -180.0 180.0 }  Krst 64.0
  GridSample :1@N :1@C8 :1@C9 :1@N1 { -180.0 180.0 }  Krst 64.0
  GridSample :1@O :1@C1 :1@N :1@C8 { -180.0 180.0 }  Krst 64.0
  GridSample :1@C1 :1@N :1@C2 :1@C3 { -180.0 180.0 }  Krst 64.0

  RandomSample :1@C2 :1@N :1@C1 {  105.0 135.0 }   Krst 256.0
  RandomSample :1@C2 :1@N :1@C8 {  105.0 135.0 }   Krst 256.0

  combine 1 2 3 4  
  count 10
  verbose 1

  % Controls on the quantum mechanical operations
  qmlev    'B3LYP',
  basis    '6-31G**',

  % Output controls: TAKE NOTE, this will generate a lot of files so
  % do this in a clean directory that you won't ls very often.
  write   'pdb', 'gaussian'
  outbase 'Conf', 'Conf'
  outsuff 'pdb',  'com'
&end
```
- create `.pdb` for structure for next step and `.com` for Gaussian input file


### 1.3.1 Sampling in `&configs`

This workflow uses two sampling strategies:

#### GridSample

Performs **uniform sampling at regular intervals** within a specified range.
Used to systematically scan dihedral angles (e.g., −180° to 180°).
- here we scan backbone Φ/Ψ/ω/χ

#### RandomSample

Performs **random sampling within a specified range**.
A flat-bottom harmonic restraint is applied, with stiffness controlled by `Krst`.

#### combine

Used only with `GridSample` operations to generate multidimensional grid sampling.
Not applied to `RandomSample`.



---
## 2. usage of `gaussian.py`
```bash
python gaussian.py
```
- align the input format for Gaussian
- it'll search every folders that include `*.com` file
- replace the setup for the first four lines

- ```bash
    1  %nprocshared=8
    2  %mem=14GB
    3  # b3lyp/6-31g(d,p)
    4
    ```
---
## 3. run calcaulation in Gaussian or ORCA
here we use Gaussian 09 in Stella
```bash
sbatch single_point.sh
```
- make sure the number of array, `#SBATCH --array=0-9`, is same as your input files
- `#SBATCH --cpus-per-task=8` and `SBATCH --mem=16G` should align your input setup 
```bash
    #!/bin/bash
    #SBATCH --job-name=g09_array
    #SBATCH --output=%x.%A_%a.out
    #SBATCH --error=%x.%A_%a.err
    #SBATCH --array=0-9
    #SBATCH --nodes=1
    #SBATCH --ntasks=1
    #SBATCH --cpus-per-task=8
    #SBATCH --mem=16G
    #SBATCH --time=08:00:00
    #SBATCH --partition=normal
```

## 4. create `.cdf` file
to make all conformers into one file via `cpptraj`

- example
```bash
parm ../../3_antechamber/S01/s01.top
trajin Conf1.pdb
trajin Conf2.pdb
trajin Conf3.pdb
trajin Conf4.pdb
trajin Conf5.pdb
trajin Conf6.pdb
trajin Conf7.pdb
trajin Conf8.pdb
trajin Conf9.pdb
trajin Conf10.pdb
trajout coords.cdf netcdf
run
quit
```
