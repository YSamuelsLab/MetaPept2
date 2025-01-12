# -*- coding: utf-8 -*-
"""
Created on Sun Apr 30 13:21:45 2023

@author: brachaer
"""


import sys
import os
import subprocess
import pandas as pd
import numpy as np
import csv


class netMHCpan_II:
	def __init__(self, args):
		# define class variables

		# define netMHCpan binary (executable) path
		self._netMHCpan_II_bin = os.path.join(args.netMHCpan_II_path, 'bin/netMHCpan_II')

		self._args = args
		self._MHC_II_list = []
		self._exp_alleles = ''

		# call class functions
		self._set_netMHCpan_II_env()
		self._get_netMHCpan_II_MHC_list() 
		self._get_experiment_HLA_II_list()

		self._netMHCpan_II_xls_file = os.path.join(self._args.output_folder, 'netMHCpan_II_binding_output.xls')

		self.affinity_df = None

		if not os.path.exists(self._netMHCpan_II_xls_file):
			self.run_netMHCpan_II()		
		else:
			print("\033[1;31m %s netMHCpan_II output file exists.." %(self._netMHCpan_II_xls_file))
			print(' to re-run netMHCpan_II delete output files from %s' %(args.output_folder))
			print('\x1b[6;30;42m' + '' + '\x1b[0m')

	def _set_netMHCpan_II_env(self):
		os.environ['TMPDIR'] = self._args.tmpdir

		# set netMHCpan_II enivronment variable
		# deafult: './tools/netMHCpanII-4.2/Linux_x86_64/'

		os.environ['NETMHCpan_II']  = self._args.netMHCpan_II_path

		if not os.path.exists(self._args.netMHCpan_II_path) or not os.path.exists(self._args.tmpdir):
			print('%s does not exist, netMHCpan_II cannot run\nwithout defining netMHCpan_II environment location' %(self._arg.netMHCpan_II_path))
			sys.exit(1) 

	def _get_netMHCpan_II_MHC_list(self):
		res = subprocess.run([self._netMHCpan_II_bin, "-listMHC_II"], capture_output=True)

		# split results into a list
		self.MHC_II_list = res.stdout.decode().splitlines()

	def _get_experiment_HLA_II_list(self):
		""" read HLA_II alleles file
		check availability in alleles file, and return list of alleles

		@args HLA_file: list of HLA alleles used in the experiment
		@type HLA_file: path to HLA text file, must be comma delimited line/s with Allele list
		example: 
		HLA-A01:01, HLA-A02:01, HLA-B4403, HLA-B5701, HLA-C0602, HLA-C1601
		"""

		# HLA_II_file = os.path.join(self._args.input_folder, self._args.alleles)
		HLA_II_file = self._args.alleles

		if not os.path.exists(HLA_II_file):
			print('%s file does not exist... exiting...' %(HLA_II_file))
			sys.exit(1)

		# read HLA_II alleles file
		# alleles should be comma delimited <- this option is not relevant, the correction is below
		# and may be in separate lines <- it does not work, the correction is below
		# self._exp_alleles = ','.join(pd.read_csv(HLA_II_file, sep="\s*,\s*", quoting=csv.QUOTE_NONE, engine='python'))
		self._exp_alleles = ','.join(pd.read_csv(HLA_II_file, header=None, engine='python')[0].values.tolist())

		# check if all given HLA alleles are available in netMHCpan
		try: 
			all(item in self._exp_alleles for item in self._MHC_II_list)
		except:
			for allele in self._exp_alleles:  
				if (not allele in self._MHC_II_list):
					print('%s not available in MHCpan_II, check  experiment HLA_II nomenclature' %(allele))
					print('find allelenames file in : ./tools/netMHCpanII-4.1\2/data/allelenames')
					print('netMHCpan_II will not work ...')
					print('exiting ...')
					sys.exit(1)

	def run_netMHCpan_II(self):
		""" 
		Run the netmhcpan command line
		Pointing directly to the binary. 
		Avoids using the tcsh script.
		
		(required setting the os.eviron manually)

		"""

		# set netMHCpan_II command line

		# [-p]                 0                    Use peptide input <netMHC_II_peptides_file>
		# peptides_file : example:  peptides.pep contains a list of peptides one peptide per line

		# [-a line]            HLA-A02:01           MHC_II allele <self._exp_alleles>
		# hla_alleles_list:  example: 'HLA-A01:01,HLA-A02:01,HLA-B4403,HLA-B5701,HLA-C0602,HLA-C1601'

		# [-BA] Include Binding affinity prediction

		netMHC_II_peptides_file =  os.path.join(self._args.output_folder, 'peptides.pep') 
		netMHCpan_II_output = open(os.path.join(self._args.output_folder, 'netMHCpan_II_binding_output.txt'), 'w')

		try:
			subprocess.run([self._netMHCpan_II_bin, "-BA", "-p", netMHC_II_peptides_file, "-a", self._exp_alleles, "-xls", "-xlsfile", self._netMHCpan_II_xls_file], stdout = netMHCpan_II_output)
			netMHCpan_II_output.close()

		except Exception as e:
			print (e)
			sys.exit(1)

		return

	# binding affinity antilog
	def _BA_exp(self,x):
		return 50000 ** (1 - x)

	def parse_netMHCpan_II(self):
		netMHCpan_II_affinity_file = os.path.join(self._args.output_folder, 'netMHCpan_II_HLA_affinity.csv')

		if os.path.exists(netMHCpan_II_affinity_file):
			print("\033[1;31m %s netMHCpan_II HLA affinity file exists.." %(netMHCpan_II_affinity_file))
			print(' to parse netMHCpan_II affinity delete output files from %s' %(self._args.output_folder))
			print('\x1b[6;30;42m' + '' + '\x1b[0m')
		else:
			# read netMHCpan_II xls output
			# collect HLA alleles data in a dictionary with binding affinity
			# [-rth float]         2             Rank Threshold for high binding peptides
			# [-rlt float]         10             Rank Threshold for low binding peptides

			# SB = < 2
			# WB = 10 - 2.0
			# NB = > 10
			
			# BA (binding affinity data) 
			# EL (eluted ligand data) 

			try:
				netMHC_II_xls_df = pd.read_csv(self._netMHCpan_II_xls_file, sep = '\t', header=[0,1])

			except (FileNotFoundError, IOError):
				print("\033[1;31m %s netMHCpan_II xls file missing.." %(self._netMHCpan_II_xls_file))
				return

			hla_cols = [ x[0] for x in netMHC_II_xls_df.columns if 'HLA_II' in str(x)]

			el_rank_df = netMHC_II_xls_df.filter(like='EL_Rank')
			el_rank_df.columns = [hla + ' EL_Rank' for hla in hla_cols]

			# find lowest EL rank
			el_rank_df.loc[:,'HLA_II rank'] = el_rank_df.values.min(axis=1)

			# find column name with lowest rank value (Lower rank == stronger binding) 
			el_rank_df.loc[:,'HLA_II Allele'] = el_rank_df.idxmin(axis=1)

			# clean the values in the new column to have only the Allele name
			el_rank_df.loc[:,'HLA_II Allele'] = el_rank_df.loc[:,'HLA_II Allele'].map(lambda value: value.split()[0])

			# set condition for ascertaining binding strength (as per netMHCpan_II default)
			m1 = el_rank_df.loc[:,'HLA_II rank'] < 2.0
			m2 = el_rank_df.loc[:,'HLA_II rank'].between(2.0, 10.0)
			el_rank_df.loc[:,'HLA_II affinity'] = np.select([m1, m2], ['SB', 'WB'], default = "")

			# add peptide column to and move to be first column
			peptide_col = netMHC_II_xls_df.loc[:, (slice(None), "Peptide")].to_numpy()
			el_rank_df.insert(loc=0, column='Peptide', value=peptide_col)

			# add affinity in nm columns

			# 50000 ** (1 - BA-score)
			aff_df = netMHC_II_xls_df.filter(like = 'BA-score')
			aff_df = aff_df.apply(self._BA_exp)
			
			aff_df.columns = [hla + ' Aff(nM)' for hla in hla_cols]

			self.affinity_df = el_rank_df.join(aff_df) 

			# save netMHC_out_df to csv
			self.affinity_df.to_csv(netMHCpan_II_affinity_file, index=False)

