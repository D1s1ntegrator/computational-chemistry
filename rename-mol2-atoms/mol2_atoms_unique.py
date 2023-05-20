#!/usr/bin/env python3
import sys

file_name = sys.argv[1]
file_name_out = sys.argv[2]

extention = ".mol2"
extention_output = ".mol2"

if not file_name.endswith(extention):
    file_name += extention

if not file_name_out.endswith(extention_output):
    file_name_out += extention_output

file = open(file_name, "r")
lines = []
for line in file.readlines():
    lines.append(line)
file.close()

# Searching for the beginning and the end of various data.
c = 0
while c < len(lines):
    if lines[c][:13] == "@<TRIPOS>ATOM":
        atoms_start = c + 1
        break
    c += 1

c += 1
while c < len(lines):
    if lines[c][:13] == "@<TRIPOS>BOND":
        atoms_end = c
        break
    c += 1

atom_names = []

# Adding data to list.
for c in range(atoms_start, atoms_end):
    line = lines[c]
    atom_names.append(line[8:12].strip())  # Location for atom names.

for c in range(len(atom_names)):
    if atom_names[c] in atom_names[:c]:
        element = atom_names[c][0]
        i = 1  # Assigning new unique name to atom.
        while i < 1000:
            new_name = element + str(i)
            if new_name not in atom_names[:c]:
                break
            i += 1
        atom_names[c] = new_name


# Replacing renamed atoms.
for c in range(atoms_start, atoms_end):
    line = lines[c]
    line = list(line)
    line[8:12] = list("%-4s" % atom_names[c - atoms_start])
    lines[c] = ''.join(line)


# Writing out .mol2 file.
file = open(file_name_out, 'w')
for line in lines:
    file.write(line)
file.close()
