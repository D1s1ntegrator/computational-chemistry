#!/usr/bin/env python3
import sys
from math import sqrt, sin, asin, degrees, radians
import numpy as np

filename1 = sys.argv[1]
filename2 = sys.argv[2]
toggle = sys.argv[2:]

if "-at" in toggle:
    doPrintAtomTypes = True
else:
    doPrintAtomTypes = False

extention = ".mol2"

if not filename1.endswith(extention):
    filename1 += extention

if not filename2.endswith(extention):
    filename2 += extention

# =============================================================================
# Functions.
# =============================================================================


def reverse_tuple(tup):
    """Reverses a tuple and returns a tuple."""
    return tuple(reversed(tup))


def add_rev_elements(ls):
    """Returns a set of tuple input elements with reversed input elements"""
    return set(ls + list(map(reverse_tuple, ls)))


def findAtomsBonds(filename):
    """Searching for atoms and bonds in .mol2 format."""
    file = open(filename, 'r')
    lines = []
    for line in file.readlines():
        lines.append(line)
    file.close()

    c = 0
    while True:
        if lines[c].startswith("@<TRIPOS>MOLECULE"):
            break
        c += 1
    c += 2
    line = lines[c].split()
    atom_count = int(line[0])
    bond_count = int(line[1])

    while True:
        if lines[c].startswith("@<TRIPOS>ATOM"):
            break
        c += 1
    c += 1
    atom_end = c + atom_count

    coordinates = {}
    names = {}
    types = {}
    while c < atom_end:
        line = lines[c].split()
        idx = int(line[0])
        names[idx] = line[1]
        x = float(line[2])
        y = float(line[3])
        z = float(line[4])
        types[idx] = line[5]
        coordinates[idx] = (x, y, z)
        c += 1

    while True:
        if lines[c].startswith("@<TRIPOS>BOND"):
            break
        c += 1
    c += 1
    bond_end = c + bond_count

    bonds = []
    while c < bond_end:
        line = lines[c].split()
        idx1 = int(line[1])
        idx2 = int(line[2])
        bonds.append((idx1, idx2))
        c += 1
    return names, coordinates, bonds, types


def calc_dist(a, b):
    """Calcualtes distance between two points."""
    dist = sqrt((a[0] - b[0])**2 + (a[1] - b[1])**2 + (a[2] - b[2])**2)
    return dist


def calc_angle(a, b, c):
    """Calcualtes angle between three points.
    Source:
    https://stackoverflow.com/questions/35176451/python-code-to-calculate-angle-between-three-point-using-their-3d-coordinates

    """
    a = np.array(a)
    b = np.array(b)
    c = np.array(c)

    ba = a - b
    bc = c - b

    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)

    return np.degrees(angle)


def calc_dih(a, b, c, d):
    """
    Calcualtes dihedal angle between four points. Uses Praxeolitic formula.
    Source:
    https://stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python

    """
    a = np.array(a)
    b = np.array(b)
    c = np.array(c)
    d = np.array(d)

    b0 = -1.0*(b - a)
    b1 = c - b
    b2 = d - c

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))


def find_connected(idx1, all_bonds):
    connected = []
    for bond in all_bonds:
        if idx1 in bond:
            idx2 = bond[1 - bond.index(idx1)]
            connected.append(idx2)
    return connected


# =============================================================================
# Determination of all existing bonds, angles and dihedrals.
# =============================================================================

atom_names, coordinates1, bonds_idx, atom_types = findAtomsBonds(filename1)
atom_names2, coordinates2, bonds_idx2, atom_types2 = findAtomsBonds(filename2)
if len(atom_names) != len(atom_names2):
    raise IOError("Atom count does not match between files!")

if atom_names != atom_names2:
    print("!!! Atom ordering does not match between files !!!")
    for a, b in zip(atom_names.values(), atom_names2.values()):
        print(a, b)
    print("Proceed with caution!")
    input()

if len(bonds_idx) != len(bonds_idx2):
    print(len(bonds_idx))
    print(len(bonds_idx2))
    raise IOError("Bond count does not match between files!")

if add_rev_elements(bonds_idx) != add_rev_elements(bonds_idx2):
    print("Bond ordering does not match between files!")
    print("Nonexisting bonds in molecule1:", set(bonds_idx2) - set(bonds_idx))
    print("Nonexisting bonds in molecule2:", set(bonds_idx) - set(bonds_idx2))
    print("Proceed with caution!")
    input()

if atom_types != atom_types2 and doPrintAtomTypes:
    print("Names of atom types does not match match between files!")
    print("First molecule atom types will be printed.")
    input()


angles_idx = []
for idx in atom_names.keys():
    connected = find_connected(idx, bonds_idx)
    for i in connected:
        for j in connected:
            if i != j:
                angle1 = (i, idx, j)
                angle2 = (j, idx, i)
                if not {angle1, angle2} & set(angles_idx):
                    angles_idx.append(angle1)


dihedrals_idx = []
for idx1, idx2 in bonds_idx:
    connected1 = find_connected(idx1, bonds_idx)
    connected1.remove(idx2)
    connected2 = find_connected(idx2, bonds_idx)
    connected2.remove(idx1)

    if connected1 == 0 or connected2 == 0:
        continue

    for i in connected1:
        for j in connected2:
            dihedral = (i, idx1, idx2, j)
            dihedrals_idx.append(dihedral)


# =============================================================================
# Calculation and printing of all existing bonds, angles and dihedrals.
# =============================================================================

S_abs = 0
S_sq = 0
print("%-7s %5s %5s %5s %5s" % ("Bond", "L1, Å", "L2, Å", "ΔL", "ΔL^2"))
for i, j in bonds_idx:
    dist1 = calc_dist(coordinates1[i], coordinates1[j])
    dist2 = calc_dist(coordinates2[i], coordinates2[j])
    bond_name = "%s-%s" % (atom_names[i], atom_names[j])

    if doPrintAtomTypes:
        bond_atom_types = "%s-%s" % (atom_types[i], atom_types[j])
    else:
        bond_atom_types = ''

    dist_diff = abs(dist1 - dist2)
    print("%-7s %5.3f %5.3f %5.3f %5.3f %s" % (
        bond_name, dist1, dist2, dist_diff, dist_diff**2, bond_atom_types
        )
    )
    S_abs += dist_diff
    S_sq += dist_diff**2

n = len(bonds_idx)
print("Total bonds:", n)
print("LAD, LAD/N, LS, LS/N")
print("%.3g, %.3g, %.3g, %.3g" % (S_abs, S_abs / n, S_sq, S_sq / n))
print()


S_abs = 0
S_sq = 0
print("%-11s %5s %5s %5s %5s" % ("Angle", "A1, °", "A2, °", "ΔA", "ΔA^2"))
for i, j, k in angles_idx:
    angl1 = calc_angle(coordinates1[i], coordinates1[j], coordinates1[k])
    angl2 = calc_angle(coordinates2[i], coordinates2[j], coordinates2[k])
    angle_name = "%s-%s-%s" % (atom_names[i], atom_names[j], atom_names[k])

    if doPrintAtomTypes:
        angl_atom_types = "%s-%s-%s" % (
            atom_types[i],
            atom_types[j],
            atom_types[k]
        )
    else:
        angl_atom_types = ''

    angl_diff = abs(angl1 - angl2)
    print("%-11s %5.1f %5.1f %5.1f %5.1f %s" % (
        angle_name, angl1, angl2, angl_diff, angl_diff**2, angl_atom_types
        )
    )
    S_abs += angl_diff
    S_sq += angl_diff**2

n = len(angles_idx)
print("Total bonds:", n)
print("LAD, LAD/N, LS, LS/N")
print("%.3g, %.3g, %.3g, %.3g" % (S_abs, S_abs / n, S_sq, S_sq / n))
print()


S_abs = 0
S_sq = 0
print("%-15s %6s %6s %6s %6s" % ("Dihedral", "D1, °", "D2, °", "ΔD", "ΔD^2"))
for i, j, k, l in dihedrals_idx:
    dih1 = calc_dih(
        coordinates1[i],
        coordinates1[j],
        coordinates1[k],
        coordinates1[l]
    )

    dih2 = calc_dih(
        coordinates2[i],
        coordinates2[j],
        coordinates2[k],
        coordinates2[l]
    )

    dihedral_name = "%s-%s-%s-%s" % (
        atom_names[i],
        atom_names[j],
        atom_names[k],
        atom_names[l]
    )
    if doPrintAtomTypes:
        dih_atom_types = "%s-%s-%s-%s" % (
            atom_types[i],
            atom_types[j],
            atom_types[k],
            atom_types[l]
        )
    else:
        dih_atom_types = ''

    dih_diff = degrees(asin(abs(sin(radians(dih1 - dih2) / 2))) * 2)
    print("%-15s %6.1f %6.1f %6.1f %7.1f %s" % (
        dihedral_name, dih1, dih2, dih_diff, dih_diff**2, dih_atom_types
        )
    )
    S_abs += dih_diff
    S_sq += dih_diff**2
n = len(dihedrals_idx)
print("Total dihedrals:", n)
print("LAD, LAD/N, LS, LS/N")
print("%.3g, %.3g, %.3g, %.3g" % (S_abs, S_abs / n, S_sq, S_sq / n))
print()
