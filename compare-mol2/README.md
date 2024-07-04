# Sybyl Mol2 molecule structure comparison

Comparison of the geometry of the same molecule is a nontrivial task. Cartesian coordinates are uninformative, and intrinsic coordinates are unambiguous. However, geometries can be fruitfully compared by considering bonds and associate lengths, angles, and dihedral angles.
This program compares molecule geometry parameters based on the bond configuration. It provides absolute deviation and squared deviation for your statistical needs. The format employed is Tripos Mol2, which you can read more about [here](http://chemyang.ccnu.edu.cn/ccb/server/AIMMS/mol2.pdf).

## Usage

The program is run from the terminal with the following arguments:
```
python mol2_atoms_unique.py input_a input_b
```

Positional arguments:<br />
`input_a` - name of the first `.mol2` input file<br />
`input_b` - name of the second `.mol2` input file<br />
