## instruciton

find the constrain atoms which is your capping group.

```bash
grep -n "constrain 0" /path/to/your/*.nw
```

use your output file to extract residue charges only into dat file.

```bash
python extract.py /path/to/your.out --show-removed -o /your/charge.dat

```
print the removed atoms

```bash
python extract.py /path/to/your.out --show-kept -o /your/charge.dat

```
print the kept atoms