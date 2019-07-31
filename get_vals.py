import sys
import numpy as np

def get_mol(filename):
  ''' read in labels and coordinates
  '''
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

def get_vec(atom1, atom2):
	return atom2 - atom1

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

def print_stuff(name, vx, vy, angle_x, angle_y, dihedral):

	print(name[:-4], vx, vy, angle_x, angle_y, dihedral)


if __name__ == "__main__":

	filename = sys.argv[1]

	numAtoms, labels, coords = get_mol(filename)

	C1 = coords[0]
	C2 = coords[1]
	X  = coords[2]
	Y  = coords[3]

	dist_x 		= get_dist(get_vec(X, C1))
	dist_y 		= get_dist(get_vec(C2, Y))

	angle_x 	= get_angle(get_vec(C1, X), get_vec(C1, C2))
	angle_y 	= get_angle(get_vec(C2, Y), get_vec(C2, C1))

	dihedral 	= get_dihedral(get_vec(X, C1), get_vec(C1, C2), get_vec(C2, Y))

	print_stuff(filename, dist_x, dist_y,  angle_x, angle_y, dihedral)
