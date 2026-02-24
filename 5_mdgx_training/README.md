# training set
- finding electron energies in every single point calculation of conformers
    - combine energies into `.dat` file
- combine all conformers into a `.cdf` file
```bash
# exampke .cdf file
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
```bash
python extract.py
```
- i'll print energies of every conformers into order