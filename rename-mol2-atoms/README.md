# Sybyl Mol2 file duplicate atom name renamer

While building a molecule using molecular building software, you can encounter a situation when you have two atoms with the same name (not to be confused with the SYBYL atom type).
Duplicate atom names create a problem in some molecule modelling software, like [Gromacs](https://www.gromacs.org/), where the force field cannot differentiate atoms in the same molecule.
For this reason, this script checks and renames the atom to a unique name.
It does this by going through all atoms sequentially: if the same atom name is detected, a new name is assigned while keeping the original element symbol and adding a unique index.
You can read more about the Tripos Mol2 file format [here](http://chemyang.ccnu.edu.cn/ccb/server/AIMMS/mol2.pdf).

## Usage

The program is run from the terminal with the following arguments:
```
python mol2_atoms_unique.py input output
```

Positional arguments:<br />
`input` - name of the `.mol2` input file<br />
`output` - name of the `.mol2` output file<br />
