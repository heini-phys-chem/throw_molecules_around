#!/usr/bin/env python3

import sys
import numpy as np
import openbabel as ob
from MDAnalysis import Universe, topology



def get_mol(filename):
	mol = ob.OBMol()
	conv = ob.OBConversion()
	conv.SetInAndOutFormats('xyz', 'xyz')
	conv.ReadFile(mol, filename)

	return mol, conv

def write_mol(mol, name):
	conv.WriteFile(mol, name)

def check_num_fragments(mols):
	if len(mols) == 2:
		return True
	else:
		return False

def check_Y_fragments(mols):
#	if mols[0].NumAtoms() == 1 and mols[0].GetFormula() in ['H', 'F', 'Cl', 'Br']:
#		return True
#	elif mols[1].NumAtoms() == 1 and mols[1].GetFormula() in ['H', 'F', 'Cl', 'Br']:
#		return True
	if mols[0].NumAtoms() > 1 and mols[1].NumAtoms() > 1:
		return False
	else:
		return True

def check_constraints_E2(mol):
	C2 = ob.OBAtom
	H  = ob.OBAtom
	Y  = ob.OBAtom

	C2 = mol.GetAtom(2)
	H  = mol.GetAtom(4)
	Y  = mol.GetAtom(mol.NumAtoms())

	angle = mol.GetAngle(C2, H, Y)

	if angle > 178 and angle < 182:
		return True
	else:
		return False

def check_constraints_SN2(mol):
	C1 = ob.OBAtom
	X  = ob.OBAtom
	Y  = ob.OBAtom

	C1 = mol.GetAtom(1)
	X  = mol.GetAtom(3)
	Y  = mol.GetAtom(mol.NumAtoms())

	angle = mol.GetAngle(X, C1, Y)

	if angle < 182 and angle > 178:
		return True
	else:
		return False

def get_non_zero(dic, key):
	if dic[key] == 0:
		return ""
	elif dic[key] == 1:
		return key
	else:
		return key + str(dic[key])

def check_sum_formula(mol, name, mols):
	target = name
	letters = target.split('_')

	if letters[5] == 'A': numElements = {'C' : 2, 'H' : 0, 'Br' : 0, 'Cl' : 0, 'F' : 0, 'N' : 0, 'O' : 0}
	if letters[5] != 'A': numElements = {'C' : 2, 'H' : 1, 'Br' : 0, 'Cl' : 0, 'F' : 0, 'N' : 0, 'O' : 0}
	#print(mol.NumAtoms())
	for i, letter in enumerate(letters):
		if i <= 3:
			if letter == 'A':
				numElements['H']  += 1
			if letter == 'B':
				numElements['N']  += 1
				numElements['O']  += 2
			if letter == 'C':
				numElements['C']  += 1
				numElements['N']  += 1
			if letter == 'D':
				numElements['C']  += 1
				numElements['H']  += 3
			if letter == 'E':
				numElements['N']  += 1
				numElements['H']  += 2
		if i == 4:
			if letter == 'A':
				numElements['F']  += 1
			if letter == 'B':
				numElements['Cl'] += 1
			if letter == 'C':
				numElements['Br'] += 1
		if i == 5:
			if letter == 'A':
				numElements['H']  += 1
			if letter == 'B':
				numElements['F']  += 1
			if letter == 'C':
				numElements['Cl'] += 1
			if letter == 'D':
				numElements['Br'] += 1

	num_formula_from_name = 'C' + str(numElements['C']) + get_non_zero(numElements, 'H') + get_non_zero(numElements, 'Br') + get_non_zero(numElements, 'Cl') + get_non_zero(numElements, 'F') + get_non_zero(numElements, 'N') + get_non_zero(numElements, 'O')

	if num_formula_from_name == mol.GetFormula():
		return True
	else:
		return False

def check_CY_dist_SN2(mol):
	dist = {1 : 1.14, 9 : 1.41, 17 : 1.86, 35 : 2.04}
	C1 = ob.OBAtom
	Y  = ob.OBAtom

	C1 = mol.GetAtom(1)
	Y  = mol.GetAtom(mol.NumAtoms())

	vec_CY = np.linalg.norm(np.array([C1.GetX() - Y.GetX(), C1.GetY() - Y.GetY(), C1.GetZ() - Y.GetZ()]))

	if vec_CY > dist[Y.GetAtomicNum()]:
		return True
	else:
		return False

def check_HY_dist_SN2(mol):
	dist = {1 : 0.78, 9 : 0.96, 17 : 1.33, 35 : 2.48}
	H  = ob.OBAtom
	Y  = ob.OBAtom
	Y  = mol.GetAtom(mol.NumAtoms())

	numAtoms = mol.NumAtoms()
	distance = 10

	for i in range(1, numAtoms):
		if mol.GetAtom(i).GetAtomicNum() != 1:
			continue
		else:
			H = mol.GetAtom(i)

			tmp_distance = np.linalg.norm([ Y.GetX() - H.GetX(), Y.GetY() - H.GetY(), Y.GetZ() - H.GetZ() ])

			if tmp_distance < distance: distance = tmp_distance

	if tmp_distance > dist[Y.GetAtomicNum()]:
		return True
	else:
		return False

def check_dist(f):
	u = Universe(f)
	print(u.atoms.positions)
	print(topology.guessers.guess_bonds(u.atoms, u.atoms.positions))
	print(len(u.atoms) - 1)

	for i in topology.guessers.guess_bonds(u.atoms, u.atoms.positions):
		if len(u.atoms) - 1 in i:
			return False
	return True

def test(mol):
	numAtoms = mol.NumAtoms()

	for i in range(1, numAtoms + 1):
		print(mol.GetAtom(i).GetAtomicNum())

if __name__ == '__main__':

	lines = open(sys.argv[1], 'r').readlines()

	for a, line in enumerate(lines):
		tokens = line.split()
		filename = tokens[0] + '.xyz'
		name = tokens[1]
		rxn = tokens[2]
		energy = float(tokens[3])

		check_dist(filename)
		exit()
		mol, conv = get_mol(filename)
		mols = mol.Separate()
		print('-------------------------------------')

		check = check_num_fragments(mols)
		if check == True: check = check_Y_fragments(mols)
		if rxn == 'e2':
			if check == True: check = check_constraints_E2(mol)
		if rxn == 'sn2':
			if check == True: check = check_CY_dist_SN2(mol)
			if check == True: check = check_HY_dist_SN2(mol)
			if check == True: check = check_constraints_SN2(mol)
			if check == True: check = check_dist(filename)
		if check == True: check = check_sum_formula(mol, name, mols)

		if check == True:
			print(filename + '\t' + name + '\t' + str(energy) + '\tok')
		else:
			print(filename + '\t' + name + '\t' + str(energy) + '\tFail')
