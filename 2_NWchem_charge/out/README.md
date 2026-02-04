# instruciton

## 🔍find the constrain atoms which is your capping group.

```bash
grep -n "constrain 0" /path/to/your/*.nw
```
command line would show:
```bash
constrain 0 #a #b #c ... #n
```
---
## 📍use your output file to extract residue charges only into dat file.

```bash
python extract.py /path/to/your.out --show-removed -o /your/charge.dat

```
command line would show:
```bash
[OK] Parsed # atoms from RESP table.
Cap atom indices to REMOVE (1-based, space/comma separated): #a #b #c ... #n

=== Removed (cap) atoms ===
   #a atom    x y z   ESP value   constr value
   #b atom    x y z   ESP value   constr value
   #c atom    x y z   ESP value   constr value
   ...
   ...
   ...
   #n atom    x y z   ESP value   constr value

=== Summary ===
Removed atoms : # of atoms
Remaining    : # of atoms
Sum(all constr)     = +0.000000
Sum(remaining constr)= value
```
code would ask you:
```bash
Round decimals? (blank = no rounding; e.g., 4):

```
press enter to leave a blank or type number wtih the decimal places you want.
command line would show:
```bash
[DONE] Wrote # charges to charge.dat
```
---
print the removed atoms

```bash
python extract.py /path/to/your.out --show-kept -o /your/charge.dat

```

print the kept atoms

---
## ✂️trim your xyz file into uncapped
```bash
python trim_xyz.py
```

* Input XYZ: `/path/to/your/xtbopt.xyz`
* Indices to remove: `1 2 3 18 19 20 15-17 27-32`
* Output: `/path/to/create/uncapped.xyz`
---
## ↔️convert xyz file into mol2 file
```bash
obabel /path/to/your.xyz -O capped.mol2
```



