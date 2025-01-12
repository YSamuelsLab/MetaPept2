#!/usr/bin/env python3

"""ContaminationSearch.py: The search of contaminated data in MSFragger outputs."""

__author__ = "Dmitry Malko"


import argparse
import re
from Bio import SeqIO
from FastFastaSearch import FastaSearch


class ContSearch:
    _sqlite_file_default = ':memory:'
    _decoy_column_name = 'Contaminant'

    def __init__(self, ref_fasta_file, msfragg_fasta_file, db_file=None):
        self._db_file = db_file if db_file else self._sqlite_file_default
        self._ref_headers = set()
        self._cont_fasta = None

        with open(ref_fasta_file) as fasta_handle:
            for header, seq in SeqIO.FastaIO.SimpleFastaParser(fasta_handle):
                self._ref_headers.add(header)

        if not re.search(r'\.fas$', msfragg_fasta_file):
            raise ValueError('ERROR: wrong input file name ({}), the file extension must be `.fas`'.format(msfragg_fasta_file))

        self._cont_fasta = re.sub(r'\.fas$', '.contaminant.fas', msfragg_fasta_file)
        with open(msfragg_fasta_file, 'r') as fasta_handle:
            with open(self._cont_fasta, 'wt') as cont_handler:
                for header, seq in SeqIO.FastaIO.SimpleFastaParser(fasta_handle):
                    if header not in self._ref_headers:
                        cont_handler.write(''.join(['>', header, '\n', seq, '\n']))

    # end of __init__()

    def filter(self, input_name, seq_column_name, keep_decoy=False):
        search_engine = FastaSearch(self._db_file)
        search_engine.set_db(self._cont_fasta, true_match=False)
        data = search_engine.db_search(input_name, seq_column_name, self._decoy_column_name, i2l_mode=False)
        if not keep_decoy:  # remove Decoy sequences from the output file
            data = data[data[self._decoy_column_name].str.len() == 0]

        return data

    # end of filter()


def main():
    parser = argparse.ArgumentParser(description='A script for the search of contaminated data in MSFragger outputs.')
    parser.add_argument('-fa', required=True, help='reference FASTA file for MSFragger')
    parser.add_argument('-fas', required=True, help='contaminated FASTA file from MSFragger')
    parser.add_argument('-i', required=True, help='CSV file with peptide sequences to remove contamination')
    parser.add_argument('-o', default='output.csv', required=False, help='output file')
    parser.add_argument('-p', default='Peptide Sequence', required=False,
                        help='query column with peptide sequences in the input CSV file')
    parser.add_argument('-db', required=False,
                        help='database file (if not specified, the data will be stored in memory')
    parser.add_argument('-k', action='store_true', help='keep `Decoy` sequences in the output file')

    args = parser.parse_args()
    fa_file = args.fa
    fas_file = args.fas
    input_file = args.i
    output_file = args.o
    seq_column_name = args.p
    db_file = args.db
    keep_decoy = args.k

    try:
        search_engine = ContSearch(fa_file, fas_file, db_file)
        data = search_engine.filter(input_file, seq_column_name, keep_decoy)
        data.to_csv(output_file, sep='\t', index=False)
    except Exception as err:
        print('Something went wrong: {}'.format(err))
    else:
        print('Contamination clean up: done')

    return None

# end of main()


if __name__ == '__main__':
    main()
