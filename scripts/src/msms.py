# -*- coding: utf-8 -*-
"""
Created on Sun Apr 16 12:44:13 2023

@author: Bracha Erlanger Avigdor, Dmitry Malko


"""


import re
import os

import pandas as pd
import numpy as np

from collections import defaultdict


class msms:
	def __init__(self, args):
		# define class variables
		self._args = args
		self._exp_design_dict = None
		self._msms_dict = defaultdict(dict)
		self._msms_df = None

		# define dataframe for each category
		self._pivot_msms_df = None
		self._msms_pivot_file = os.path.join(self._args.output_folder, 'pivot_msms.csv')

		if not os.path.exists(self._msms_pivot_file):
			print("processing msms data")
			# call class functions
			self._read_experimental_design()
			self._parse_msms()
			self._pivot_msms_df.to_csv(self._msms_pivot_file, index=False)
		else:
			print("\033[1;31m %s pivot msms output file exists.." %(self._msms_pivot_file))
			print(' to re-run msms module delete output files from %s' %(args.output_folder))
			print('\x1b[6;30;42m' + '' + '\x1b[0m')

			self._pivot_msms_df = pd.read_csv(self._msms_pivot_file, sep=',', engine='python')

	def _calc_fragmentation(self, fragmentation_ions ,Sequence):
		""" calculates spectra fragmentations

		@args fragmentation_ions: a string of ions from msms spectra
		@type fragmentation_ions: 
		example: 
		'y1;y2;y3;y4;y5;y6;y8;a2;b2;b3;b4;b5;b6;b7;b8;b9'

		returns a list of the breaks

		"""
		# remove background ions, i.e. y1-H2O, y1-NH3
		ions = [x for x in fragmentation_ions.split(';') if '-' not in x]

		# extract y and b ions
		ions = re.findall(r'[yb]\d+', ','.join(ions))

		# and collapse redundant ions
		ions = np.unique(ions).tolist()

		breaks = [0] * len(Sequence)

		# indexes = [1 for i, x in enumerate(ions) if x.find('y')]

		# iterate over ions list
		for i,x in enumerate(ions):
			# check for carboxyl terminus (C-term) ions
			# get the position of the break and reverse for y break points
			if 'y' in x:
				try:
					y_index = - int(x.split('y')[1])
					breaks[y_index] += 1
				except:
					pass

			# check for amino terminus (N-term) ions
			# get the position of the break and keep original position
			if 'b' in x:
				try:
					b_index = int(x.split('b')[1])
					breaks[b_index] += 1
				except:
					continue

		return breaks

	# calculate coverage
	def _calc_coverage(self, breaks, Sequence):
		""" calculates %fragmentation coverage

		@args breaks: list of spectra fragmentation breaks
		@type

		@args Sequence: amino acid sequence
		@type

		"""

		# count # of fragmentation breaks
		breaks_count = np.count_nonzero(breaks)

		# calculate coveage as:
		# The ratio of fragmentation breaks to possible breaks in the sequence
		# possible breaks == length of sequence - 1

		coverage = round(100 * (breaks_count/(len(Sequence) - 1)), 2)

		return coverage

	def calc_consecutive_no_breaks(self):

		return

	def _parse_msms(self):
		msms_file = os.path.join(self._args.input_folder, self._args.msms_file)
		self._msms_df = pd.read_csv(msms_file, sep='\t', engine='python')
		self._msms_df = self._msms_df[self._msms_df['Length'] >= 8]  # HARD CODDED FILTER!!!

		# ions = msms_df['Matches'][0]
		# Sequence = msms_df['Sequence'][0]

		# breaks = calc_fragmentation(ions, Sequence)
		# coverage = calc_coverage(breaks, Sequence)
		# load msms file and select only desired columns
		colnames = [
			'Raw file',
			'Sequence',
			'Scan number',
			'Length',
			'Charge',
			'Mass',
			'Retention time',
			'PEP',
			'Score',
			'Delta score',
			'Matches',
			'Intensities'
		]

		if 'Expectation' in self._msms_df.columns:  # for MSFragger data
			colnames = [
				'Raw file',
				'Sequence',
				'Scan number',
				'Length',
				'Charge',
				'Mass',
				'Retention time',
				'PeptideProphet',
				'Hyperscore',
				'Delta score',
				'Expectation'
			]


		self._msms_df = self._msms_df.loc[:, colnames]

		# change raw file name with experiment name from  experimental design file
		self._msms_df.replace({'Raw file': self._exp_design_dict}, inplace = True)

		if 'Hyperscore' not in self._msms_df.columns:
			# MSFragger has no data to calculate coverage - this branch only for MaxQuant

			# calculate fragmentation breaks for each peptide in the dataframe
			self._msms_df['breaks'] = self._msms_df.apply(lambda row: self._calc_fragmentation(row['Matches'], row['Sequence']), axis=1)

			# calculate coverage for each peptide in the dataframe
			self._msms_df['coverage'] = self._msms_df.apply(lambda row: self._calc_coverage(row['breaks'], row['Sequence']), axis=1)

		self._msms_df_get_max_score()

		self._msms_df_get_Scan_numbers()
		return

	def _msms_df_get_max_score(self):
		# group by sequence and experiment
		# get sample with best MaxQuant Score
		if 'Score' in self._msms_df.columns:
			# MaxQuant data
			self._msms_df = self._msms_df.loc[self._msms_df.groupby(['Sequence', 'Raw file'])['Score'].idxmax()]
		elif 'Hyperscore' in self._msms_df.columns:
			# MSFragger data
			self._msms_df = self._msms_df.loc[self._msms_df.groupby(['Sequence', 'Raw file'])['Hyperscore'].idxmax()]

		return

	def _make_column_names_from_exp_data(self, col_text):
		col_names = [v + ' ' + col_text for v in self._exp_design_dict.values()]

	def _msms_df_get_Scan_numbers(self):
		colnames = [
			'Scan number',
			'Score',
			'Delta score',
			'PEP',
			'Charge',
			'Mass',
			'Retention time',
			'Matches',
			'Intensities',
			'breaks',
			'coverage'
		]

		if 'Expectation' in self._msms_df.columns:  # for MSFragger data
			colnames = [
				'Scan number',
				'Hyperscore',
				'Delta score',
				'PeptideProphet',
				'Charge',
				'Mass',
				'Retention time',
				'Expectation'
			]

		self._pivot_msms_df = pd.pivot_table(self._msms_df,
				values=colnames,
				index='Sequence',
				columns='Raw file',
				sort=False).reset_index()

		self._pivot_msms_df.columns = [' '.join(col).strip() for col in self._pivot_msms_df.columns.values]

	def _read_experimental_design(self):
		# exp_design_file = os.path.join(self._args.input_folder, self._args.exp_design_file)
		exp_design_file = self._args.sample_desc_file
		exp_design_df = pd.read_csv(exp_design_file, sep='\t', engine='python')
		exp_design_df['Source_File'] = exp_design_df['Source_File'].apply(lambda x: re.sub(r'\.[^.]+$', '', x))

		# convert the experimental design dataframe to dictionary
		# key: experiment raw file name
		# value: experiment description
		self._exp_design_dict = dict(zip(exp_design_df.Source_File, exp_design_df.Experiment))

