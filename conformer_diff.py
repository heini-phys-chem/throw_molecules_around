#!/usr/bin/env python3

import sys
import numpy as np
import operator

def read_file(f):
	lines = open(f, 'r').readlines()

	d = dict()

	for line in lines:
		tokens = line.split()
		tmp = tokens[0].split('/')
		name = tmp[1][:-1] + tmp[3].split('-')[3]
		energy = float(tokens[2])*627.509

		if name in d:
			d[name].append([energy])
		else:
			d[name] = [[energy]]

	return d

def get_max_diff(d):

	diffs = dict()

	for key, value in d.items():
		diffs[key] = np.abs(np.max(value) - np.min(value))

	return diffs


if __name__ == "__main__":
	filename = sys.argv[1]

	d = read_file(filename)
	d = get_max_diff(d)

	print(max(d.items(), key=operator.itemgetter(1)))
