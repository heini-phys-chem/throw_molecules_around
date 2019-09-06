#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style('whitegrid')
sns.set_style('ticks')
sns.set_context("poster")

fs = 25
ls = 20
legend = 18
ms = 55
lw = .4

def read_data(filename):
	return np.load(filename)

def get_data(filenames, rxn):
	if rxn == 'X':
		C_X = read_data(filenames[0][:-1])
		F_X  = read_data(filenames[1][:-1])
		Cl_X = read_data(filenames[2][:-1])
		Br_X = read_data(filenames[3][:-1])

		return [C_X, F_X, Cl_X, Br_X]

	elif rxn == 'Y':
		C_Y = read_data(filenames[0][:-1])
		H    = read_data(filenames[1][:-1])
		H_Y  = read_data(filenames[2][:-1])
		F_Y  = read_data(filenames[3][:-1])
		Cl_Y = read_data(filenames[4][:-1])
		Br_Y = read_data(filenames[5][:-1])

		return [C_Y, H, H_Y, F_Y, Cl_Y, Br_Y]

	else:
		C_Y = read_data(filenames[0][:-1])
		H_Y  = read_data(filenames[1][:-1])
		F_Y  = read_data(filenames[2][:-1])
		Cl_Y = read_data(filenames[3][:-1])
		Br_Y = read_data(filenames[4][:-1])

		return [C_Y, H_Y, F_Y, Cl_Y, Br_Y]



def make_figure():
	fig = plt.figure(figsize=(10,8))

	ax1 = plt.subplot(212)
	ax2 = plt.subplot(221)
	ax3 = plt.subplot(222)

	axes = [ax1, ax2, ax3]

	fig.text(0.5, 0.03, r'Distance [$\mathrm{\AA}$]', ha='center', va='center', fontsize=fs)
	fig.text(0.04, 0.5, r'Distance [$\mathrm{\AA}$]', ha='center', va='center', rotation='vertical', fontsize=fs)

	for i in range(len(axes)):
		axes[i].tick_params(axis='both', labelsize=ls)

	return fig, axes

def scatter(data, ax, colors, labels, marker, rxn):

	for i, atom in enumerate(data):
		if rxn == "E2":
			x = [ j[0] for j in atom]
			y = [ j[1] for j in atom]
			if i == 0:
				sns.scatterplot(x,y, ax=ax, color=colors[i], marker=marker[i], s=ms, linewidth=0)
			else:
				sns.scatterplot(x,y, ax=ax, color=colors[i], marker=marker[i], s=ms, linewidth=lw)
		elif rxn == "SN2":
			atom = [i for i in atom if i[0] > -2]
			x = [ j[0] for j in atom]
			y = [ j[1] for j in atom]
			if i == 0:
				sns.scatterplot(y,x, ax=ax, color=colors[i], marker=marker[i], s=ms, linewidth=0)
			else:
				sns.scatterplot(y,x, ax=ax, color=colors[i], marker=marker[i], s=ms, linewidth=lw)

	plt.legend(labels=labels, bbox_to_anchor=(-0.1, 1.30), loc='upper center', ncol=6, fontsize=legend, fancybox=True, shadow=True)


if __name__ == "__main__":

	filenames_X     = open(sys.argv[1], 'r').readlines()
	filenames_Y     = open(sys.argv[2], 'r').readlines()
	filenames_sn2   = open(sys.argv[3], 'r').readlines()
	filenames_sn2_Y = open(sys.argv[4], 'r').readlines()

	X        = get_data(filenames_X, 'X')
	Y        = get_data(filenames_Y, 'Y')
	sn2      = get_data(filenames_sn2, 'X')
	sn2_Y    = get_data(filenames_sn2_Y, 'SN2')

	colors_X = ['k', 'C0', 'C1', 'C2']
	labels_X = ['C', 'F', 'Cl', 'Br']
	marker_X = ['o', 'v', '>', 'd']

	colors_Y = ['k', 'gray', 'C3', 'C0', 'C1', 'C2']
	labels_Y = ['C', 'H', r'H$^-$', 'F$^-$', 'Cl$^-$', 'Br$^-$']
	marker_Y = ['o', '^', 's' ,'v', '>', 'd']

	c_sn2_Y = ['k', 'C3', 'C0', 'C1', 'C2']
	l_sn2_Y = ['C', r'H$^-$', 'F$^-$', 'Cl$^-$', 'Br$^-$']
	m_sn2_Y = ['o', 's' ,'v', '>', 'd']

	fig, axes = make_figure()

	scatter(X,     axes[1], colors_X, labels_X, marker_X, 'E2')
	scatter(Y,     axes[2], colors_Y, labels_Y, marker_Y, 'E2')
	scatter(sn2,   axes[0], colors_X, labels_X, marker_X, 'SN2')
	scatter(sn2_Y, axes[0], c_sn2_Y, l_sn2_Y, m_sn2_Y, 'SN2')

	fig.savefig('TS_geometries.png')
	fig.savefig('TS_geometries.pdf')

	#plt.show()
