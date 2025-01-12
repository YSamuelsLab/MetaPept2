# -*- coding: utf-8 -*-
"""
Created on Wed June 07 12:56:53 2023

@author: Bracha Erlanger Avigdor, Dmitry Malko

"""


import pandas as pd
import numpy as np
import re
import os
import sys
sys.path.insert(0, './src')


# Class of different styles
class style():
    BLACK = '\033[30m'
    RED = '\033[31m'
    GREEN = '\033[32m'
    YELLOW = '\033[33m'
    BLUE = '\033[34m'
    MAGENTA = '\033[35m'
    CYAN = '\033[36m'
    WHITE = '\033[37m'
    UNDERLINE = '\033[4m'
    RESET = '\033[0m'


class filterTables:
    def __init__(self, args):
        # define class variables
        self._args = args

    def _print_file_exists_message(self, filename):
        print("%s %s file exists.." %(style.BLUE, filename))
        print(' to re-run the module, delete from output folder:  %s' %(self._args.output_folder))
        print(style.RESET)

    def _print_file_not_exists(self, filename, message=""):
        print("%s %s file does not exists.." %(style.RED, filename))
        print(message)
        print(style.RESET)

    def _set_peptide_priority(self, df):
        # we assign a negtive number for "hits" according to condition
        # where the permutation with highest priority is assigned with the lowest negative number.
        # we will then select the sequence with the lowest "hit" number to filter out permutations that
        # do not take precedence over the original sequence.

        conditions = [
            (df['Permutation_Index']==0) & (df['CDS'].notna()),  # -5
            (df['Permutation_Index']==0) & ((df['CDS'].isna()) & (df['nuORFs'].notna())),  # -3
            (df['Permutation_Index']==0) & (df[['CDS', 'nuORFs']].isna().all(1)),  # -1
            (df['Permutation_Index'] > 0) & (df['CDS'].notna()),  # -4
            (df['Permutation_Index'] > 0) & ((df['CDS'].isna()) & (df['nuORFs'].notna()))  # -2
            ]

        choices = [-5, -3, -1, -4, -2]

        return(np.select(conditions, choices, default = df['Permutation_Index']))

    def filter_MQ_netMHCpan_peptides(self):
        # check if files exist:
        IMP_unfiltered_file = os.path.join(self._args.output_folder, 'IMP_unfiltered.csv')

        if not os.path.exists(IMP_unfiltered_file):
            self._print_file_not_exists(IMP_unfiltered_file, "%s is essential for IMP" %(IMP_unfiltered_file))
            return

        IMP_filtered_file = os.path.join(self._args.output_folder, 'IMP_filtered.csv')

        if os.path.exists(IMP_filtered_file):
            self._print_file_exists_message(IMP_filtered_file)
            return

        # read unfiltered IMP file
        df = pd.read_csv(IMP_unfiltered_file, sep=',', engine='python')

        # define database 'hits' status
        df['hits'] = self._set_peptide_priority(df)

        # perform some extra filtering steps
        df = (df[df['hits'] == df.groupby('Sequence')['hits'].transform(min)]
            .drop(df.filter(regex='[a-zA-Z\s]Count').columns, axis=1)
            .drop(df.filter(regex=re.compile(r'amino acid', re.IGNORECASE)).columns, axis=1)
            .query('`HLA affinity` in ["SB", "WB"]')
            .query('Length >= 8 & Length <= {}'.format(self._args.max_len)))

        # Remove the columns from their original positions.
        # Add the columns back but in desired positions.
        df.insert(1, 'Sequence_Permutations', df.pop('Sequence_Permutations'))
        df.insert(2, 'Permutation_Index', df.pop('Permutation_Index'))
        df.insert(3, 'CDS', df.pop('CDS'))
        df.insert(4, 'nuORFs', df.pop('nuORFs'))
        df.insert(5, 'hits', df.pop('hits'))

        df.to_csv(IMP_filtered_file, index=False)

