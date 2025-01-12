# -*- coding: utf-8 -*-
"""
Created on Wed June 07 12:56:53 2023

@author: Bracha Erlanger Avigdor

"""


import pandas as pd
import os
import sys


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


class mergeTables:
    def __init__(self, args):
        # define class variables
        self._args = args

    def _print_file_exists_message(self, filename):
        print("%s %s file exists.." %(style.BLUE, filename))
        print(' to re-run the module, delete from output folder:  %s' %(self._args.output_folder))
        print(style.RESET)

    def _print_file_not_exists(self, filename, message = ""):
        print("%s %s file does not exists.." %(style.RED, filename))
        print(message)
        print(style.RESET)

    def merge_MQ_netMHCpan_tables(self):
        # combine peptide, msms and netMHCpan output
        IMP_unfiltered_file = os.path.join(self._args.output_folder, 'IMP_unfiltered.csv')
        # read msms_pivot table

        if os.path.exists(IMP_unfiltered_file):
            self._print_file_exists_message(IMP_unfiltered_file)
            return

        msms_pivot_file = os.path.join(self._args.output_folder, 'pivot_msms.csv')
        msms_pivot_df = pd.read_csv(msms_pivot_file, sep=',',engine='python')

        # read netMHCpan affinity table
        HLA_aff_file = os.path.join(self._args.output_folder, 'netMHCpan_HLA_affinity.csv')
        HLA_affinity_df = pd.read_csv(HLA_aff_file, sep=',',engine='python')

        IL_peptides_file = os.path.join(self._args.output_folder, 'IL_peptides.csv')
        peptides_df = pd.read_csv(IL_peptides_file, sep=',',engine='python')

        IMP_unfiltered_df = peptides_df.merge(HLA_affinity_df, left_on='Sequence_Permutations', right_on='Peptide', how='left').merge(msms_pivot_df, on='Sequence', how='left')
        # df = pd.merge(self._peptides_df, HLA_affinity_df,left_on='Sequence_Permutations', right_on='Peptide', how='left')
        # IMP_unfiltered_df = pd.merge(df, msms_pivot_df, on='Sequence', how='left')

        IMP_unfiltered_df.to_csv(IMP_unfiltered_file, index=False)

