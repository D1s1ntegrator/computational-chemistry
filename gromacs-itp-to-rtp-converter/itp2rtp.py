#!/usr/bin/env python3
import sys

file_name = sys.argv[1]

extention = ".itp"
extention_output = ".rtp"

if not file_name.endswith(extention):
    file_name += extention

file_name_out = file_name[:-len(extention)] + extention_output

# To delete atoms, specify file with their names in the terminal.
if len(sys.argv) > 2:
    file = open(sys.argv[2], "r")
    lines = []
    for line in file.readlines():
        lines.append(line.split())
        lines = list(filter(None, lines))
    atoms2delete = set(sum(lines, []))
    file.close()
    print("Deleted atoms:", ", ".join(atoms2delete))
else:
    atoms2delete = set()

file = open(file_name, "r")
lines = []
for line in file.readlines():
    lines.append(line.split())
lines = list(filter(None, lines))
file.close()

# Searching for the beginning and the end of various data.
for i, line in enumerate(lines):
    if line == ['[', 'atoms', ']']:
        atoms_start = i
    if line == ['[', 'bonds', ']']:
        bonds_start = i
    if line == ['[', 'pairs', ']']:
        pairs_start = i
        break

atoms_start += 3
atoms_end = bonds_start
bonds_start += 2
bonds_end = pairs_start

atom_names_full = []
atom_names = []
atom_types = []
charge = []
charge_group = []
precision = 0

# Adding data to appropriate lists.
for c in range(atoms_start, atoms_end):
    line = lines[c]
    atom_names_full.append(line[4])
    if not line[4] in atoms2delete:
        atom_names.append(line[4])
        atom_types.append(line[1])
        charge.append(float(line[6]))
        precision = max(len(line[6].split(".")[1]), precision)
        charge_group.append(int(line[5]))


# Reordering charge groups after deleting atoms.
new_charge_group = [1, ]
for i in range(1, len(charge_group)):
    if charge_group[i] == charge_group[i - 1]:
        new_charge_group.append(new_charge_group[-1])
    else:
        new_charge_group.append(new_charge_group[-1] + 1)
charge_group = new_charge_group


bond_i = []
bond_j = []
# Finding bonded atom names by indexes.
for c in range(bonds_start, bonds_end):
    line = lines[c]
    i = int(line[0]) - 1
    j = int(line[1]) - 1
    if not bool({atom_names_full[i], atom_names_full[j]} & atoms2delete):
        bond_i.append(atom_names_full[i])
        bond_j.append(atom_names_full[j])

# Writing out .rtp file.
file = open(file_name_out, 'w')
file.write("[ " + lines[3][0] + " ]\n")
file.write("  [ atoms ]\n")

for out in zip(atom_names, atom_types, charge, charge_group):
    # file.write("{:>9s}{:>9s}{:>8.4f}{:>4d}\n".format(*out))
    file.write("%9s%9s%8.4f%4d\n" % out)

file.write("  [ bonds ]\n")

for out in zip(bond_i, bond_j):
    file.write("%9s%6s\n" % out)

charge = round(sum(charge), precision * 2) + 0
print(("Total charge: %." + str(precision) + "f") % charge)
print("Total atoms: %d" % len(atom_names))
print("Total bonds: %d" % len(bond_i))
file.write("\n")
file.close()
