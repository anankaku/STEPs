# STEPsv1.0
Force field parameters, build tools, coordinate and energy files, and sample input files for Yet Another Peptoid force field (STEPs).

# Requirements to use parameter files:
* ambertools, tested with ambertools22
* Python to use sampleParmed.py to convert to Gromacs
* An MD engine, i.e: Gromacs, Amber, OpenMM etc.

# Requirements to follow parameterization steps:
* OpenBabel, produced with version 3.1.0
* xTB, produced with version 6.5.1
* ORCA, produced with ORCA 5.0.3 but tested on ORCA 4.0
* NWChem
* ambertools, tested with ambertools22

# Installation
Install prerequisites as needed. Unzip STEPs.zip and place in desired location.

# Usage (STEPs)
Refit parameters are stored in `BuildTools/Parameters/STEPsv1.dat`
Additional OTB cap group parameters are stored in `BuildTools/Parameters/addOTB.dat`
tleap input file `BuildTools/leaprc.q4mdfft` contains the four (4) cap groups and seventy (70) residues. 
To use tleap input file, modify line 467 of `BuildTools/leaprc.q4mdfft` (optionally modify line 465 by commenting out the OTB parameters) to build peptoid from sequence.
Modify lines 469 and 470 with desired output names
execute `tleap -f leaprc.q4mdfft` to produce desired output files 

# Additional tleap Information
More in depth details for building structures from sequence can be found in the AmberTools manual.
Specific bond, angles, and dihedrals can be imposed on the structure prior to the saveamberparm line as described in the AmberTools manual.
Common Peptoid Dihedrals and examples
Omega: `impose peptoid { 2 } {{ "C1" "N" "C" "O" 0.00 }}` sets capped monomer to cis conformation
Phi: `impose peptoid { 1 2 3 } {{ "C" "N" "CA" "C" 180.00 }}` sets backbone phi dihedral to 180 degrees
Psi: `impose peptoid { 1 2 3 } {{ "N" "CA" "C" "N" 180.00 }}` sets backbone psi dihedral to 180 degrees
Chi1: `impose peptoid { 2 } {{ "C" "CA" "N" "C1" 180.00 }}` sets capped monomer chi1 dihedral to 180 degrees

# Usage (Parameter Fitting)
More specific details for the mdgx fitting procedure can be found in the amber specific tutorial here https://ambermd.org/tutorials/advanced/tutorial32/index.php

To generate new parameters a sample mdgx input file is provided in `FittingExample/fitting.mdgx`. This file may depend on the slightly modified gaff2.dat file which is included in this location (a few dihedrals had to be added to the default gaff2.dat to appease mdgx, they were refit anyways and are included in the resulting output parameter file, making it primarily a concern only at this step).

Coordinates and energies for two iterations as described in the manuscript are included in the folder `Coords/` Each residue additionally has a pdb version of the resulting structures from the first round of mdgx structure generation stored in the `Coords/RXX/rama/` folders. These coords.cdf and energies.dat files are the ones used in the mdgx fitting procedure

The included sample mdgx input file can be run with `mdgx -i fitting.mdgx` 

# Usage (Generation of additional side chains)
* In depth details of the procedure are described in the manuscript and and the tutorial linked above

As an example of the process for R20 (N-benzyl) several of the sample files needed can be found in the included `SampleInputs` folder.

Briefly:
* Beginning in the folder `SampleInputs/0_StructureOpt/`
The capped structure for R20 can be generated in OpenBabel using the command `obabel -:"CC(=O)N(Cc1ccccc1)CC(=O)N(C)C" -ismi -oxyz -O R20.xyz -h --gen3d`
This resulting `R20.xyz` can then be used as the xtb input (input file `SampleInputs/0_StructureOpt/R20/xtb.orca`) and executed with `orca xtb.orca > R20_xtb.out` in that folder
* For NWChem RESP Charge generation `SampleInputs/1_RESP/`
The final coordinates from `SampleInputs/0_StructureOpt/xtb.orca.xyz` are then used in sample NWChem input file `SampleInputs/1_RESP/R20/R20.nw` sample constraints can be seen at the bottom of the file. Outputs are included in `SampleInputs/1_RESP/R20/outputs`
* For mdgx orca structure generation:
Convert `xtb.orca.xyz` to pdb using vmd -> resulting file `RXXvmd.pdb`. 
- `antechamber -i RXXvmd.pdb -fi pdb -c rc -cf R20charge.dat -o R20.prepi -fo prepi`
- Use `tleap -f GenTopology.tleap` to generate `R20.top` and `R20.crd`
- use RXXvmd for referencing scan coordinates in the orca_structures file
- generate input structures for orca `SampleInputs/2_ambermdgx/R20/orca_structures/gen_coor.mdgx`
- `mdgx -i gen_coor.mdgx`
- use Orca to calculate energies
- grep energies and coordinates and use for fitting procedure outlined above
