# training set
- finding electron energies in every single point calculation of conformers
    - combine energies into `.dat` file
- combine all conformers into a `.cdf` file
    - check `.cdf` frames
    - `cpptraj -p ../3_antechamber/S01/s01.top -y S01/coords.cdf -tl`
```bash
# example .cdf file
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

## find electron energy of every conformers
to make energies into one file
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