import sys
import numpy as np


def print_test_dist(coords, place, d_C1Xs, d_C2Ys):
	if place == True:
		print("X dist_old:\t\t%.5f" % np.linalg.norm(get_vec(coords[0], coords[2])))
		print("Y Dist_old:\t\t%.5f" % np.linalg.norm(get_vec(coords[1], coords[3])))
		print()
	else:
		print("X dist should be:\t%.5f" % d_C1Xs)
		print("X dist_new:\t\t%.5f" % np.linalg.norm(get_vec(coords[0], coords[2])))
		print("Y Dist should be:\t%.5f" % d_C2Ys)
		print("Y Dist_new:\t\t%.5f" % np.linalg.norm(get_vec(coords[1], coords[3])))
		print()


def print_test_angle(coords, place, a_C1Xs, a_C2Ys):
	if place == True:

		print("X angle_old:\t\t%.5f" % get_angle(get_vec(coords[0], coords[2]),
				                                     get_vec(coords[0], coords[1])))
		print("Y angle_old:\t\t%.5f" % get_angle(get_vec(coords[1], coords[3]),
				                                     get_vec(coords[1], coords[0])))

		print()
	else:
		print("X how it should be:\t%.5f" % a_C1Xs)
		print("X angle_new:\t\t%.5f" % get_angle(get_vec(coords[0], coords[1]),
				                                     get_vec(coords[0], coords[2])))
		print("Y how it should be:\t%.5f" % a_C2Ys)
		print("Y angle_new:\t\t%.5f" % get_angle(get_vec(coords[1], coords[3]),
					                                   get_vec(coords[1], coords[0])))
		print()

def print_test_dihedral(coords, place, dihedrals):
	if place == True:

		print("dihedral_old:\t\t%.5f" % get_dihedral(
				get_vec(coords[2], coords[0]), 
				get_vec(coords[0], coords[1]),
				get_vec(coords[1], coords[3]) ))
	else:
		print("how it should be:\t%.5f" % dihedrals)
		print("dihedral_new:\t\t%.5f" % get_dihedral(
				get_vec(coords[2], coords[0]), 
				get_vec(coords[0], coords[1]),
				get_vec(coords[1], coords[3]) ))

def get_data(filename):
	''' read in new distances, angles and dihedral'''
	lines = open(filename, 'r').readlines()

	fins, fouts, d_C1Xs, d_C2Ys, a_C1Xs, a_C2Ys, dihedrals = [ [] for i in range(7) ]

	for line in lines:
		tokens = line.split()
		fin = "xyz/" + tokens[0] + '.xyz'
		fout = "xyz_out/" +tokens[0] + '_2.xyz'
		d_C1X = float(tokens[1])
		d_C2Y = float(tokens[2])
		a_C1X = float(tokens[3])
		a_C2Y = float(tokens[4])
		dihedral = float(tokens[5])

		fins.append(fin)
		fouts.append(fout)

		d_C1Xs.append(d_C1X)
		d_C2Ys.append(d_C2Y)
		a_C1Xs.append(a_C1X)
		a_C2Ys.append(a_C2Y)
		dihedrals.append(dihedral)

	return fins, fouts, d_C1Xs, d_C2Ys, a_C1Xs, a_C2Ys, dihedrals

def get_mol(filename):
	''' read in labels and coordinates'''
	lines = open(filename, 'r').readlines()

	labels = []
	coords = []

	for i, line in enumerate(lines):
		if i == 0: numAtoms = line
		if i == 1: continue

		if i > 1:
			tokens = line.split()
			labels.append(tokens[0])
			coords.append(np.array([ float(tokens[1]), float(tokens[2]), float(tokens[3]) ]))

	return numAtoms, labels, coords

def center_mol(atom, coords):
  ''' center mol '''

  for i in range(len(coords)):
    coords[i] = [ coords[i][0] - atom[0], coords[i][1] - atom[1], coords[i][2] - atom[2] ]

  return coords

def get_vec(atom1, atom2):
	''' return a vector from 2 points '''
	return np.array(atom2) - np.array(atom1)

def get_angle(v1, v2):
	cos = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
	arccos = np.arccos(cos)
	degrees = np.degrees(arccos)

	return degrees

def get_dihedral(v1, v2, v3):
	v12 = np.cross(v1, v2)
	v23 = np.cross(v2, v3)

	n12 = v12 / np.linalg.norm(v12)
	n23 = v23 / np.linalg.norm(v23)

	torsion = get_angle(n12, n23)

	return torsion

def get_dist(v):
	return np.linalg.norm(v)

def rotation_matrix(axis, theta):
	''' rotational matrix '''
	axis = np.asarray(axis)
	axis = axis / np.sqrt(np.dot(axis, axis))
	a = np.cos(theta / 2.0)

	b, c, d = axis * np.sin(theta / 2.0)
	aa, bb, cc, dd = a * a, b * b, c * c, d * d
	bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d

	return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
									 [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
									 [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])


def print_xyz(numAtoms, labels, coords):
	print(str(numAtoms))
	for i, coord in enumerate(coords):
		print("%s\t%.4f\t%.4f\t%.4f" % (labels[i], coord[0], coord[1], coord[2]))

def move_atom_along_vector(origin, atom, v, d_new, coords, XY):
	d_old = get_dist(v)
	vn = v / np.linalg.norm(v)
	diff = d_new - d_old
	if XY == 'X': coords[2] = atom + diff * vn
	if XY == 'Y': coords[3] = atom + diff * vn

	return coords 

def rotate_along_angle(angle_new, v1, v2, coords, XY):
	v12 = np.cross(v1, v2)
	n12 = v12 / np.linalg.norm(v12)

	angle_old = get_angle(v1, v2)

	axis = n12
	theta = np.radians(angle_new - angle_old)

	if XY == 'X':
		coords = center_mol(coords[0], coords)
		atom = coords[2]
		coords[2] = np.dot( rotation_matrix(axis, theta), [atom[0], atom[1], atom[2]] )
	if XY == 'Y':
		coords = center_mol(coords[1], coords)
		atom = coords[3]
		coords[3] = np.dot( rotation_matrix(axis, theta), [atom[0], atom[1], atom[2]] )

	return coords

def rotate_along_dihedral(atom, dihedral_new, v1, v2, v3, coords):
	dihedral_old = get_dihedral(v1, v2, v3)
	diff = dihedral_old - dihedral_new

	axis = v2
	theta = np.radians(diff)

	coords = center_mol(coords[1], coords)

	coords[3] = np.dot( rotation_matrix(axis, theta), [atom[0], atom[1], atom[2]] )

	return coords

def get_atoms(coords):
	C1 = coords[0]
	C2 = coords[1]
	X  = coords[2]
	Y  = coords[3]

	return C1, C2, X, Y

def write_xyz(numAtoms, labels, corods, fout):
	f = open(fout, 'a')
	f.write(str(numAtoms) + '\n')
	for i in range(len(labels)):
		f.write("{0}\t{1:.5f}\t{2:.5f}\t{3:.5f}\n".format(labels[i], coords[i][0], coords[i][1], coords[i][2]) )
	f.close()


if __name__ == '__main__':

	filename = sys.argv[1]

	''' read in input xyz files, output file names, new distances, angles and dihedrals '''
	fins, fouts, d_C1Xs, d_C2Ys, a_C1Xs, a_C2Ys, dihedrals = get_data(filename) 

	for i, fin in enumerate(fins):
		''' loop over all molecules '''
		numAtoms, labels, coords = get_mol(fin)

		C1, C2, X, Y = get_atoms(coords)

		''' set new distances '''
		print_test_dist(coords, True, d_C1Xs[i], d_C2Ys[i])
		coords = move_atom_along_vector(C1, X, get_vec(C1, X), d_C1Xs[i], coords, 'X')
		C1, C2, X, Y = get_atoms(coords)
		coords = move_atom_along_vector(C2, Y, get_vec(C2, Y), d_C2Ys[i], coords, 'Y')
		C1, C2, X, Y = get_atoms(coords)
		print_test_dist(coords, False, d_C1Xs[i], d_C2Ys[i])

		''' set new angles '''
		print_test_angle(coords, True, a_C1Xs[i], a_C2Ys[i])
		coords = rotate_along_angle(a_C1Xs[i], get_vec(C1, X), get_vec(C1, C2), coords, 'X')
		C1, C2, X, Y = get_atoms(coords)
		coords = rotate_along_angle(a_C2Ys[i], get_vec(C2, Y), get_vec(C2, C1), coords, 'Y')
		C1, C2, X, Y = get_atoms(coords)
		print_test_angle(coords, False, a_C1Xs[i], a_C2Ys[i])


		''' set new dihedral '''
		print_test_dihedral(coords, True, dihedrals[i])
		coords = rotate_along_dihedral(Y,
				                              dihedrals[i],
																			get_vec(X, C1),
																			get_vec(C1, C2),
																			get_vec(C2, Y),
																			coords)
		print_test_dihedral(coords, False, dihedrals[i])

		''' print  new coordinates '''
#		print_xyz(numAtoms, labels, coords)
		print('\n#####################################################\n')
		write_xyz(numAtoms, labels, coords, fouts[i])
