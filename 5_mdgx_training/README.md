# training set
- finding electron energies in every single point calculation of conformers
    - combine energies into `.dat` file
- combine all conformers into a `.cdf` file
    - check `.cdf` frames
        - `cpptraj -p ../3_antechamber/S01/s01.top -y S01/coords.cdf -tl`
---
use `pre_cdf.py` to combine all conformers into a `.cdf` file
- change file path in `pre_cdf.py`
```bash
# Path to the topology file (used in cpptraj as "parm")
top_path = Path("../../3_antechamber/S0X/s0X.top")

# Path where make_cdf.in will be written
output_file = Path("make_cdf.in")

# Directory containing Conf*.pdb files
pdb_dir = Path("../../4_mdgx_topology/S0X")
```
- and follow the instruction in terminal after run `pre_cdf.py`

```bash
# example .in file
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

---
## find electron energy of every conformers
use `extract_energy.py` to make energies into one file
- `python extract_energy.py your/energy/log/folder`
- example print out
```bash
[OK] Conf1.log: -826.123654599000 Hartree
[OK] Conf2.log: -826.123654599000 Hartree
[OK] Conf3.log: -826.123654599000 Hartree
[OK] Conf4.log: -826.123654598000 Hartree
[OK] Conf5.log: -826.123654599000 Hartree
[OK] Conf6.log: -826.123654599000 Hartree
[OK] Conf7.log: -826.123654598000 Hartree
[OK] Conf8.log: -826.123654599000 Hartree
[OK] Conf9.log: -826.123654599000 Hartree
[OK] Conf10.log: -826.123654599000 Hartree

Wrote 10 energies -> /home/tuu61186/STEPs/5_mdgx_training/S01/energy.dat
```
- example output
```bash
-826.123654599000
-826.123654599000
-826.123654599000
-826.123654598000
-826.123654599000
-826.123654599000
-826.123654598000
-826.123654599000
-826.123654599000
-826.123654599000
```