# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 11:56:53 2023

@author: Bracha Erlanger Avigdor, Dmitry Malko

"""

import os
import sys
import numpy as np
import pandas as pd
import re
import itertools

# this package enables using R packages in python
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri as rpyn
import rpy2.rinterface_lib.callbacks
rpy2.rinterface_lib.callbacks.consolewrite_warnerror = lambda *args: None  # to suppress warnings from rpy2

sys.path.append('..')
from FastFastaSearch import FastaSearch


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


class peptides:
    def __init__(self, args):
        # define class variables
        self._args = args
        # define netMHCpan binary (executable) path
        self._peptides_df = None
        self.IMP_unfiltered_df = None

        # make protViz R package instance
        # for adding ssrc hydrophobicity column
        self._protviz = importr('protViz')

        # call class functions
        # read files to dataframes:

        print("reading peptides file")

        # read MQ peptides file 'peptides.txt'
        self._read_peptides()

        # read IL_peptides.csv. if not exist, call self._I_to_L
        self._process_IL_peptides()

        # write 'peptides.pep' input for netMHCpan
        self._write_petides_for_netMHCpan()

    def _print_file_exists_message(self, filename):
        print("%s %s file exists.." % (style.BLUE, filename))
        print(' to re-run the module, delete from output folder:  %s' % self._args.output_folder)
        print(style.RESET)

    def _print_file_not_exists(self, filename, message=""):
        print("%s %s file does not exists.." %(style.RED, filename))
        print(message)
        print(style.RESET)

    def _add_fasta_db_search_column(self, column_name, db_name, db_fasta_file):
        """ adds a column marking the existence of a peptide in a given fasta database

        @args column_name: Name of column containing the peptides
        @type column_name: str

        @args db_name: Name of database to be searched
        @type db_name: str

        @args db_fasta_file: File name containing the fasta database to be searched
        @type db_fasta_file: str

        """

        # make fasta database instance with class fastaDB from DB_search module
        print('performing %s database search' %(db_name))
        db_fasta_file = os.path.join(self._args.database_folder, db_fasta_file)

        """ THE OLD SLOW CODE FOR FULL-TEXT SEARCH
        fasta_db = fasta_DB(db_name, db_fasta_file)

        # apply sequence search on each peptide in collumn <column_name> of dataframe <df>
        # if exists in the database, make as true in a new column  <db_name>

        # calculate fragmentation breaks for each peptide in the dataframe
        self._peptides_df[str(db_name)] = self._peptides_df[str(column_name)].apply(fasta_db.fasta_DB_search)
        """

        search_engine = FastaSearch()
        search_engine.set_db(db_fasta_file, true_match=False)
        self._peptides_df = search_engine.db_search(self._peptides_df, str(column_name), str(db_name), i2l_mode=False)

    def _add_ssrc_column(self, column_name):
        """ adds a column with ssrc measure for hydrophobicity
        based on the method described in (Krokhin, Craig, Spicer, Ens, Standing, Beavis, and Wilkins 2004)
        from the R package 'protViz'

        protViz must be installed in the environment through R.
        here we import an R package using the rpy2 python package

        @args column_name: Name of column containing the peptides
        @type column_name: str

        """

        # apply sequence search on each peptide in column <column_name> of dataframe <df>
        # if exists in the database, make as true in a new column  <db_name>

        # calculate fragmentation breaks for each peptide in the dataframe

        # THE OLD CODE BELOW CASES A PROBLEM WITH COMPUTING EFFICIENCY!!!
        # prtotViz function should not be run on each peptide
        # self._peptides_df['ssrc_hydrophobicity'] = self._peptides_df[str(column_name)].apply(self._protviz.ssrc).str[0]

        pept_list = list(self._peptides_df[str(column_name)])
        data_ssrc = pd.DataFrame(
            {str(column_name): pept_list, 'ssrc_hydrophobicity': rpyn.rpy2py(self._protviz.ssrc(pept_list))},
            index=None)
        self._peptides_df = pd.merge(self._peptides_df, data_ssrc, on=[str(column_name)], how='left')

    def _read_peptides(self):
        peptides_file = os.path.join(self._args.input_folder, self._args.peptides_file)
        if os.path.exists(peptides_file):
            self._peptides_df = pd.read_csv(peptides_file, sep='\t', engine='python')
            self._peptides_df = self._peptides_df[self._peptides_df['Length'] >= 8]  # HARD CODDED FILTER!!!

            # Filter out peptides dataframe for Reverse/Decoy peptides
            if 'Reverse' in self._peptides_df.columns:
                # MaxQuant data
                self._peptides_df = self._peptides_df.loc[(self._peptides_df['Reverse'].astype(str) != '+') &
                    (self._peptides_df['Leading razor protein'].str.startswith('REV') is False) &
                    (self._peptides_df['Potential contaminant'].astype(str) != '+'), :]
            elif 'Contaminant' in self._peptides_df.columns:
                # MSFragger data
                self._peptides_df = self._peptides_df[self._peptides_df['Contaminant'].isna()]
            else:
                raise ValueError('ERROR: can not determine the data type: MaxQuant or MSFragger')
        else:
            self._print_file_not_exists(peptides_file, "%s is essential for IMP" %(self._args.peptides_file))
            sys.exit(1)

    def _process_IL_peptides(self):
        self._IL_peptides_file = os.path.join(self._args.output_folder, 'IL_peptides.csv')

        if os.path.exists(self._IL_peptides_file):

            self._print_file_exists_message(self._IL_peptides_file)
            self._peptides_df = pd.read_csv(self._IL_peptides_file, engine='python')
        else:

            self._print_file_not_exists(self._IL_peptides_file)
            self._I_to_L()

            # perform database search on all peptides including I to L

            # add ssrc hydrophobicity column
            self._add_ssrc_column('Sequence_Permutations')

            # add nuORFdb database search
            self._add_fasta_db_search_column('Sequence_Permutations', 'nuORFs', self._args.nuORFdb_fasta_file)

            # add CDS database search
            self._add_fasta_db_search_column('Sequence_Permutations', 'CDS', self._args.CDS_fasta_file)

            # write peptides_df to cvs file
            self._peptides_df.to_csv(self._IL_peptides_file, index=False)

    def _I_to_L_permutations(self, Sequence):
        """ Generates I to L permutations from a given amino acid sequence

        @args Sequence: amino acid sequence
        @type Sequence: String

        """

        # for each peptide in peptides_df, search for I and L
        # create a new identical line, except the I to L conversion

        # find all instances of I or L
        res = re.finditer(r'[IL]', Sequence)
        # retrieve index of I or L occurrences
        IL_index = [m.start() for m in res]

        if len(IL_index) < 1:
            return(None)

        IL_permutations = []

        # generate all possible permutations of I and L
        # in the length of occurrences of I and/or L in the peptide
        for row in itertools.product(['I', 'L'], repeat=len(IL_index)):
            s_lst = np.array(tuple(Sequence))
            # s_lst[[IL_index]] = list(row)  # FutureWarning: Using a non-tuple sequence for multidimensional indexing is deprecated
            s_lst[tuple([IL_index])] = list(row)
            IL_permutations.append(''.join(s_lst))

        # remove the original permutation
        # to ensure that the original sequence is in the first position
        IL_permutations.remove(Sequence)

        all_IL_perm = [Sequence] + IL_permutations

        return all_IL_perm

    def _I_to_L(self):
        """ manages I to L permutations in a given peptides dataframe

        @args peptides_df: peptides df derived from peptides.txt
        @type peptides_df: pandas dataframe

        """

        # add permutation index column
        # original sequence is marked 0
        # consequent I to L permutations are marked sequentially 1-n

        self._peptides_df['Sequence_Permutations'] = self._peptides_df['Sequence']
        self._peptides_df['Permutation_Index'] = 0

        dfs = []  # define list of dataframes

        for i in range(len(self._peptides_df.index)):
            row = self._peptides_df.iloc[i]
            Sequence = row['Sequence']

            IL_permutations = self._I_to_L_permutations(Sequence)

            if IL_permutations != None:
                permutation_count = len(IL_permutations)
                rep_df = pd.concat([row.to_frame().T] * permutation_count)
                rep_df['Sequence_Permutations'] = IL_permutations
                rep_df['Permutation_Index'] = list(range(permutation_count))

                # append permutation dataframe to dfs
                # remove first row, which is the original row (rep_df.tail(-1)).

                dfs.append(rep_df.tail(-1))

        # combine all permutation dataframes
        combined_dfs = pd.concat(dfs, ignore_index= True)

        # concatenate with original peptides dataframe
        self._peptides_df = pd.concat([self._peptides_df, combined_dfs], ignore_index=True)

        # sort peptides dataframe by (original) Sequence and permutation index
        self._peptides_df.sort_values(['Sequence', 'Permutation_Index'], inplace=True, ignore_index=True)

    def _write_petides_for_netMHCpan(self):
        print('writing peptides.pep for netMHCpan')
        netMHCpan_peptides_file = os.path.join(self._args.output_folder, 'peptides.pep')
        self._peptides_df['Sequence_Permutations'].to_csv(netMHCpan_peptides_file, header=False, index=False)

    def filter_peptides(self):
        # check if files exist:
        IMP_unfiltered_file = os.path.join(self._args.output_folder, 'IMP_unfiltered.csv')

        if not os.path.exists(IMP_unfiltered_file):
            self._print_file_not_exists(IMP_unfiltered_file, "%s is essential for IMP" % IMP_unfiltered_file)
            return

        IMP_filtered_file = os.path.join(self._args.output_folder, 'IMP_filtered.csv')

        if os.path.exists(IMP_filtered_file):
            self._print_file_exists_message(IMP_filtered_file)
            return

        # read unfiltered IMP file
        df = pd.read_csv(IMP_unfiltered_file, sep=',', engine='python')

        # define database 'hits' status
        df['hits'] = self._get_I_to_L_hits(df)

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
