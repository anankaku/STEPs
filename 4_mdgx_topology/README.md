# create dihedral conformers
## print mol2 2D structure figure to assign the atoms for dihedral

```bash
python mol2_atomname.py -i /path/to/s01.mol2 -o s01_names.png --hide-h-labels
```
- create png file without hydrogen
```bash
python mol2_atomname.py -i s01.mol2 -o s01_names.svg
```
- or svg file