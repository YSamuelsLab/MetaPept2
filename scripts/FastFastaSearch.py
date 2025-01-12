#!/usr/bin/env python3

"""FastFastaSearch.py: The fast search of short sequences in FASTA files."""

__author__ = "Dmitry Malko"


import argparse
import re
import os
import sys
import csv
import pandas as pd
import numpy as np
import sqlite3
from Bio import SeqIO
from CSVtools import CSV


class FastaSearch:
    _IL_dict = {'I': 'L', 'L': 'I'}

    _sqlite_file_default = ':memory:'

    _query_name = 'pept'
    _query_column_seq = 'Sequence'
    _query_column_pep = 'Query'
    _query_columns = [_query_column_seq, _query_column_pep]

    _sbjct_name = 'intersec'
    _sbjct_column_seq = 'Sequence'
    _sbjct_column_pep = 'Query'
    _sbjct_column_id = 'IDs'
    _sbjct_columns = [_sbjct_column_seq, _sbjct_column_pep, _sbjct_column_id]

    _db_name = 'db'
    _db_name_fts5 = 'db_fts5'
    _db_column_id = 'ID'
    _db_column_seq = 'Sequence'
    _db_column_header = 'Header'
    _db_columns = [_db_column_id, _db_column_seq, _db_column_header]

    def __init__(self, db_file=None):
        db_file = db_file if db_file else self._sqlite_file_default
        csv.field_size_limit(sys.maxsize)  # IMPORTANT! to open a CSV file with long fields

        if db_file != self._sqlite_file_default and os.path.exists(db_file):
            # os.remove(db_file)
            os.rename(db_file, db_file + '.bak')

        self._conn = sqlite3.connect(db_file)
        self._true_match = None

    # end of __init__()

    def set_db(self, fasta_file_name, true_match=False):
        fasta_data = []
        self._true_match = true_match
        with open(fasta_file_name) as fasta_handle:
            for header, seq in SeqIO.FastaIO.SimpleFastaParser(fasta_handle):
                seq_id = re.sub(r'^(\S+).*', r'\1', header)
                fasta_data.append((seq_id, seq, header))

        fasta_data = pd.DataFrame(fasta_data, columns=self._db_columns)
        fasta_data.to_sql(name=self._db_name, con=self._conn)

        self._conn.execute("""CREATE INDEX id_index ON {}({});""".format(self._db_name, self._db_column_id))

        if self._true_match:
            self._conn.execute("""CREATE INDEX seq_index ON {}({});""".format(self._db_name, self._db_column_seq))
        else:
            # let's use SQLite FTS5 extension to make full-text search faster (USE TRIGRAM TOKEN!)
            self._conn.execute(
                """CREATE VIRTUAL TABLE {} USING fts5 ({}, tokenize="trigram");""".format(
                    self._db_name_fts5, ','.join(self._db_columns)))
            self._conn.execute(  # copy data to fts5 table (trigram tokens don't work with content='table')
                """INSERT INTO {} ({}) SELECT {} FROM {};""".format(
                    self._db_name_fts5, ','.join(self._db_columns), ','.join(self._db_columns), self._db_name))

        self._conn.commit()

        return fasta_data

    # end of set_db()

    def i2l(self, pep):
        all_perm = set()  # a list to hold all current peptide permutations initialized with original peptide
        all_perm.add(pep)

        aa = np.array(list(pep))  # break string to list of amino acids
        il_index = np.isin(aa, ['I', 'L']).nonzero()[0]  # find index of I\L
        for i in il_index:  # for each I or L aa
            i_perm = []  # list of this i permutations
            for perm in all_perm:  # iterate over existing permutations
                # insert the current replacement of i aa, prefix of permuted pep and suffix of original pep
                new_pep = perm[:i] + self._IL_dict[aa[i]] + pep[i + 1:]
                i_perm.append(new_pep)
            all_perm.update(i_perm)  # append all i permutation to peptide permutations

        # all_perm.remove(pep)
        # pep_list = [pd.Series([pep, x], index=self._query_columns) for x in all_perm]
        # return pep_list

        return all_perm

    # end of i2l()

    def db_search(self, data, query_column, new_column, i2l_mode=False):
        if isinstance(data, str):  # if the input is a file name
            data = pd.read_csv(data, sep=CSV.get_delimiter(data), engine='python')

        if query_column not in data.columns:
            raise ValueError('ERROR: there is no column {} in the input data'.format(query_column))

        data2search = data[[query_column]].rename(columns={query_column: self._query_column_seq})
        if i2l_mode:
            data2search[self._query_column_pep] = data2search[self._query_column_seq].apply(
                lambda x: self.i2l(x)
            )
            data2search = data2search.explode(self._query_column_pep, ignore_index=False)
        else:
            data2search[self._query_column_pep] = data2search[self._query_column_seq]

        data2search.to_sql(name=self._query_name, con=self._conn, index=False)

        if self._true_match is None:
            raise ValueError('ERROR: set database before using db_search() method!')

        if self._true_match:
            query = """SELECT {}, {}, (SELECT GROUP_CONCAT(DISTINCT {} || ":" || {}) 
                FROM {} WHERE {} = {}) AS id 
                FROM {} WHERE id IS NOT NULL;""".format(
                '.'.join([self._query_name, self._query_column_seq]),
                '.'.join([self._query_name, self._query_column_pep]),

                '.'.join([self._db_name, self._db_column_seq]),
                '.'.join([self._db_name, self._db_column_id]),

                self._db_name,
                '.'.join([self._db_name, self._db_column_seq]),
                '.'.join([self._query_name, self._query_column_pep]),
                self._query_name
            )
        else:
            query = """SELECT {}, {}, (SELECT GROUP_CONCAT(DISTINCT
                    {} || "_" || INSTR({}, {}) || "/" || LENGTH({})  || ":" || {}) 
                    FROM {} INNER JOIN {} USING ({})
                    WHERE {} MATCH "{}:" || {}) AS id 
                    FROM {} WHERE id IS NOT NULL;""".format(
                '.'.join([self._query_name, self._query_column_seq]),
                '.'.join([self._query_name, self._query_column_pep]),

                '.'.join([self._query_name, self._query_column_pep]),  # PEPTIDE_SEQ
                '.'.join([self._db_name, self._db_column_seq]),  # INSTR
                '.'.join([self._query_name, self._query_column_pep]),  # INSTR
                '.'.join([self._db_name, self._db_column_seq]),  # LENGTH
                '.'.join([self._db_name, self._db_column_id]),  # ID

                self._db_name_fts5,  # FROM
                self._db_name,  # INNER JOIN
                self._db_column_id,  # USING

                self._db_name_fts5,  # WHERE
                self._db_column_seq,  # MATCH - to search only in the specified column (fts5 default: in all columns)
                '.'.join([self._query_name, self._query_column_pep]),

                self._query_name
            )

        q_data = self._conn.execute(query).fetchall()
        q_data = pd.DataFrame(q_data, columns=self._sbjct_columns)
        q_data = q_data[[self._sbjct_column_seq, self._sbjct_column_id]].groupby(self._sbjct_column_seq).agg(';'.join).reset_index()
        q_data = q_data.rename(columns={self._query_column_seq: query_column, self._sbjct_column_id: new_column})
        data = data.merge(q_data, on=[query_column], how='left')
        data[new_column] = data[new_column].fillna('')

        return data

    # end of db_search()

# end of class FastaSearch


def main():
    parser = argparse.ArgumentParser(description='A script for peptide searching in FASTA files')
    parser.add_argument('-i', required=True, help='input CSV file with query sequences')
    parser.add_argument('-p', default='Sequence', required=False,
                        help='query column with peptide sequences in the input CSV file')
    parser.add_argument('-n', required=True, help='a new column name to write the results of searching')
    parser.add_argument('-f', required=True, help='the FASTA file to search in it')
    parser.add_argument('-m', action='store_true',
                        help='direct match between peptides (default: search protein substrings')
    parser.add_argument('-o', default='output.csv', required=False, help='output file')
    parser.add_argument('-i2l', action='store_true', help='search with the I2L replacement')
    parser.add_argument('-db', required=False,
                        help='database file (if not specified, the data will be stored in memory')

    args = parser.parse_args()
    input_file = args.i
    seq_column_name = args.p
    new_column_name = args.n
    fasta_file = args.f
    true_match = args.m
    output_file = args.o
    i2l_mode = args.i2l
    db_file = args.db

    try:
        search_engine = FastaSearch(db_file)
        print('creating database...')
        search_engine.set_db(fasta_file, true_match)
        print('searching...')
        data = search_engine.db_search(input_file, seq_column_name, new_column_name, i2l_mode)
        print('saving results...')
        data.to_csv(output_file, sep='\t', index=False)
    except Exception as err:
        print('Something went wrong: {}'.format(err))
    else:
        print('Search in FASTA: done')

    return None

# end of main()


if __name__ == '__main__':
    main()
