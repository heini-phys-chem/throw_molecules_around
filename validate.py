#!/usr/bin/env python3

import sys
import numpy as np
import openbabel as ob



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
	if mols[0].NumAtoms() == 1 and mols[0].GetFormula() in ['H', 'F', 'Cl', 'Br']:
		return True
	elif mols[1].NumAtoms() == 1 and mols[1].GetFormula() in ['H', 'F', 'Cl', 'Br']:
		return True
	else:
		return False

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

	#print('##### ' + num_formula_from_name)
	#print('##### ' + mol.GetFormula())
	if num_formula_from_name == mol.GetFormula():
		return True
	else:
		return False

if __name__ == '__main__':

	lines = open(sys.argv[1], 'r').readlines()

	for a, line in enumerate(lines):
		tokens = line.split()
		filename = tokens[0] + '.xyz'
		name = tokens[1]
		rxn = tokens[2]

		mol, conv = get_mol(filename)
		mols = mol.Separate()

		check = True
		if check == True: check = check_num_fragments(mols)
#		if check == False: print('###### 1')
		if check == True: check = check_Y_fragments(mols)
#		if check == False: print('###### 2')
		if check == True:
			if rxn == 'e2':  check = check_constraints_E2(mol)
			if rxn == 'sn2': check = check_constraints_SN2(mol)
#		if check == False: print('###### 3')
		if check == True:  check = check_sum_formula(mol, name, mols)
#		if check == False: print('###### 4')

		if check == True:
			print(filename + '\t' + name + '\tok')
		else:
			print(filename + '\t' + name + '\tFail')
