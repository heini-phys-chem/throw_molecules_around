import sys
import random
import numpy as np
import matplotlib.pyplot as plt

import matplotlib.patches as mpatches

fs = 15
ls = 15

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
			coords.append([ float(tokens[1]), float(tokens[2]), float(tokens[3]) ])

	return numAtoms, labels, coords

def center_C1(coords):
	''' move C1 to origin
	'''
	C1 = coords[0]

	for i in range(len(coords)):
		coords[i] = [ coords[i][0] - C1[0], coords[i][1] - C1[1], coords[i][2] - C1[2] ]

	return coords

def rotate_C2(coords):
	''' rotate C-C bond to x-axis
	'''
	C1 = coords[0]
	C2 = coords[1]
	ex = np.array([1.0, 0.0, 0.0])

	vec = np.array([ C2[0] - C1[0], C2[1] - C1[1], C2[2] - C1[2] ])

	theta =	np.arccos( np.dot(vec, ex) / (np.linalg.norm(vec) * np.linalg.norm(ex)) )

	axis = np.cross(vec, ex)

	for i in range(len(coords)):
		coords[i] =	np.dot( rotation_matrix(axis, theta), [ coords[i][0], coords[i][1], coords[i][2] ] )

	return coords


def rotate_X(coords):
	''' rotate X into xz plane
	'''
	C1 = coords[0]
	C2 = coords[1]
	X	= coords[2]

	n1 = np.array([0,-1,0])

	vecC1X	= np.array([ X[0] - C1[0], X[1] - C1[1], X[2] - C1[2]	])
	vecC1C2 = np.array([ C2[0] - C1[0], C2[1] - C1[1], C2[2] - C1[2] ])

	n2 = np.cross( vecC1X, vecC1C2 )

	ex = np.array([1.0, 0.0, 0.0])

	if X[1] < 0:	theta =	np.arccos( np.dot(n1, n2) / (np.linalg.norm(n1) * np.linalg.norm(n2)) )
	if X[1] > 0:	theta =	-np.arccos( np.dot(n1, n2) / (np.linalg.norm(n1) * np.linalg.norm(n2)) )

	axis = ex

	for i in range(len(coords)):
		coords[i] =	np.dot( rotation_matrix(axis, theta), [ coords[i][0], coords[i][1], coords[i][2] ] )

	return coords

def center_C2(coords):
	''' move C1 to origin
	'''
	C2 = coords[1]

	for i in range(len(coords)):
		coords[i] = [ coords[i][0] - C2[0], coords[i][1] - C2[1], coords[i][2] - C2[2] ]

	return coords

def rotate_C1(coords):
	''' rotate C-C bond to x-axis
	'''
	C1 = coords[0]
	C2 = coords[1]
	ex = np.array([1.0, 0.0, 0.0])

	vec = np.array([ C2[0] - C1[0], C2[1] - C1[1], C2[2] - C1[2] ])

	theta =	np.arccos( np.dot(vec, ex) / (np.linalg.norm(vec) * np.linalg.norm(ex)) )

	axis = np.cross(vec, ex)

	for i in range(len(coords)):
		coords[i] =	np.dot( rotation_matrix(axis, theta), [ coords[i][0], coords[i][1], coords[i][2] ] )

	return coords


def rotate_Y(coords):
	''' rotate X into xz plane
	'''
	C1 = coords[0]
	C2 = coords[1]
	Y	= coords[3]

	n1 = np.array([0,1,0])

	vecC1Y	= np.array([ Y[0] - C1[0], Y[1] - C1[1], Y[2] - C1[2]	])
	vecC1C2 = np.array([ C2[0] - C1[0], C2[1] - C1[1], C2[2] - C1[2] ])

	n2 = np.cross( vecC1Y, vecC1C2 )

	ex = np.array([1.0, 0.0, 0.0])

	if Y[1] < 0:	theta =	-np.arccos( np.dot(n1, n2) / (np.linalg.norm(n1) * np.linalg.norm(n2)) )
	if Y[1] > 0:	theta =	np.arccos( np.dot(n1, n2) / (np.linalg.norm(n1) * np.linalg.norm(n2)) )

	axis = ex

	for i in range(len(coords)):
		coords[i] =	np.dot( rotation_matrix(axis, theta), [ coords[i][0], coords[i][1], coords[i][2] ] )

	return coords


def rotation_matrix(axis, theta):
	''' rotational matrix
	'''
	axis = np.asarray(axis)
	axis = axis / np.sqrt(np.dot(axis, axis))
	a = np.cos(theta / 2.0)

	b, c, d = -axis * np.sin(theta / 2.0)
	aa, bb, cc, dd = a * a, b * b, c * c, d * d
	bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d

	return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
									 [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
									 [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

def print_xyz(numAtoms, labels, coords):
	''' print xyz file
	'''
	print(numAtoms)

	for i in range(len(coords)):
		print(labels[i], coords[i][0], coords[i][1], coords[i][2])

def print_mol(coords, labels, numAtoms):
	print(numAtoms)

	for i in range(len(labels)):
		print(labels[i], coords[i][0], coords[i][1], coords[i][2])

def merge(x, y):
	lst = []

	for i in range(len(x)):
		lst.append((x[i], y[i]))

	return lst

def store(data, name):
	np.save(name, data)

if __name__ == "__main__":
	filenames =	open(sys.argv[1], 'r').readlines()
	c1x, c1y, c2x, c2y, XFx, XFy, XClx, XCly, XBrx, XBry, Hx, Hy, Hmx, Hmy = [ [] for i in range(14) ]

	for i, filename in enumerate(filenames):

		numAtoms, labels, coords = get_mol(filename[:-1])

		center_C1(coords)
		rotate_C2(coords)
		rotate_X(coords)
		#coords = center_C2(coords)
		#coords = rotate_C1(coords)
		#coords = rotate_Y(coords)

		C1 = coords[0]
		C2 = coords[1]
		X	= coords[2]
		Y	= coords[3]

		print_mol(coords, labels, numAtoms)

		c1x.append(C1[0])
		c1y.append(C1[2])
		c2x.append(C2[0])
		c2y.append(C2[2])
#		Hx.append(coords[4][0])
#		Hy.append(coords[4][2])

#		if filename[-6] == "A":
#			Hmx.append(Y[0])
#			Hmy.append(Y[2])
		if filename[-8] == "A":
			XFx.append(X[0])
			XFy.append(X[2])
		if filename[-8] == "B":
			XClx.append(X[0])
			XCly.append(X[2])
		if filename[-8] == "C":
			XBrx.append(X[0])
			XBry.append(X[2])

#	plot_the_stuff(c1x, c1y, c2x, c2y, XFx, XFy, XClx, XCly, XBrx, XBry)#, Hmx, Hmy, Hx, Hy)
	C1 = merge(c1x, c1y)
	C2 = merge(c2x, c2y)
	C = C1 + C2
#	H  = merge(Hx, Hy)
#	Hm = merge(Hmx, Hmy)
	F	 = merge(XFx, XFy)
	Cl = merge(XClx,XCly)
	Br = merge(XBrx,XBry)

	store(C,   'npy/C_X.npy')
	store(F,	 'npy/F_X.npy')
	store(Cl, 'npy/Cl_X.npy')
	store(Br, 'npy/Br_X.npy')
#	store(H,   'npy/H_Y.npy')
#	store(Hm, 'npy/Hm_Y.npy')
