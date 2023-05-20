# Gromacs itp to rtp file format converter

Protein modelling with the ligand requires a specific force field for that particular ligand.
For this reason, there exist general force files like [CGenFF](https://cgenff.umaryland.edu/). In this example, this tool generates `.str` file format, which for [Gromacs](https://www.gromacs.org/) it is necessary to convert to [Include topology](https://manual.gromacs.org/documentation/current/reference-manual/file-formats.html#itp)
`.itp` file format using [`cgenff_charmm2gmx_py3_nx2.py`](https://d1s1ntegrator.github.io/cgenff_charmm2gmx_py3_nx2.py)
script provided by MacKerell lab [webpage](https://www.charmm.org/archive/charmm/resources/charmm-force-fields/).
However, if the user wants to add new residue to the force field itself `.itp` has to be converted to [Residue topology](https://manual.gromacs.org/documentation/current/reference-manual/file-formats.html#rtp)
`.rtp` file format, where this script comes into play.
If you want to know more about adding new residue to Gromacs, you can read [here](https://www.researchgate.net/publication/323226378_Bonded_Force_Constant_Derivation_of_Lysine-Arginine_Cross-linked_Advanced_Glycation_End-Products)
or [here](https://manual.gromacs.org/current/how-to/topology.html).

## Usage

The program is run from the terminal with the following arguments:
```
python itp2rtp.py input
```

Positional arguments:<br />
`input` - name of the `.itp` input file<br />

After finishing it creates an `.rtp` file with the same filename.
