#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style('whitegrid', {'axes.spines.bottom':True})
sns.set_style('ticks')
sns.set_context("poster")


def get_df(f, ref):
	df = pd.read_csv(f, delim_whitespace=True)

	if ref == "old":
		df.columns = ["Reaction", "Calculation", "Target", "E_ts", "Rank_ts", "conf_id", "dist", "lvl", "E_react_old", "rank_react", "barrier"]
		df = df[["Target", "E_ts", "Rank_ts","E_react_old", "rank_react"]]

		df['conf_rank'] = df.groupby('Target')['E_react_old'].rank(ascending=0)
		df = df.loc[df['conf_rank'] == 1.0]

		return df[["Target", "E_ts", "Rank_ts","E_react_old"]]

	if ref == "new":
		df.columns = ["xyz", "Target", "E_react_new", "Validation"]
		df = df[["Target", "E_react_new", "xyz"]]

		df['conf_rank'] = df.groupby('Target')['E_react_new'].rank(ascending=0)
		df = df.loc[df['conf_rank'] == 1.0]

		return df[["Target", "E_react_new", "xyz"]]

def histogram(data, label):
	sns.distplot(data, kde=False, label=label)

def make_figure():
	fig, axes = plt.subplots(1, 1)
	fig.set_size_inches(8, 6, forward=True)
	plt.gcf().subplots_adjust(bottom=.20, left=.1, wspace=.25)

	return fig, axes

def set_stuff():
	fig.text(.5, .04, r"Energy [kcal$\cdot$mol$^{-1}$]", ha='center', va='center', fontsize=25)

	plt.yscale('log')
	plt.legend(fontsize=25, frameon=False)


if __name__ == "__main__":
	f_old = sys.argv[1]
	f_new = sys.argv[2]

	df_old = get_df(f_old, 'old')
	df_new = get_df(f_new, 'new')

	df = pd.merge(df_old, df_new, on='Target', how='inner')

	df['Ea_old'] = (df['E_ts'] - df['E_react_old'])*627.509
	df['Ea_new'] = (df['E_ts'] - df['E_react_new'])*627.509

	print(df)
	th_upper = 100
	th_lower = 0

	df_old = df.query('Ea_old < @th_upper & Ea_old > @th_lower')
	df_new = df.query('Ea_new < @th_upper & Ea_new > @th_lower')

	df_old.to_csv('df_old.csv', index=None, sep=' ', header=False)
	df_new.to_csv('df_new.csv', index=None, sep=' ', header=False)

	fig, axes = make_figure()

	labels = [r"S$_N$2 old", r"S$_N$2 new"]

	histogram(df_old['Ea_old'], labels[0])
	histogram(df_new['Ea_new'], labels[1])

	set_stuff()

	plt.show()
