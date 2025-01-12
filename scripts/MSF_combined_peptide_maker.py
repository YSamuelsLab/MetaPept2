#!/usr/bin/env python3

"""MSF_combined_peptide_maker.py: A script to create combined_peptide.tsv file based on psm.tsv file
if you do not have one (MSFragger does not report the file for `one replica` run)"""

__author__ = "Dmitry Malko"


import re
import argparse
from glob import glob
import pandas as pd
from CSVtools import CSV


columns = [
    'Peptide',
    'Protein Start',
    'Protein End',
    'Prev AA',
    'Next AA',
    'Peptide Length',
    'Charge',
    'Intensity',
    'Protein',
    'Protein ID',
    'Entry Name',
    'Gene',
    'Protein Description',
    'Mapped Genes',
    'Mapped Proteins',
    'Hyperscore'
]


def main():
    parser = argparse.ArgumentParser(
        description='A script to create combined_peptide.tsv file based on psm.tsv MSFragger file'
    )
    parser.add_argument('-i', required=True, help='the input psm.tsv file')
    parser.add_argument('-o', required=True, help='the output combined_peptide.tsv file')
    parser.add_argument('-s', required=True, help='the sample description file')

    args = parser.parse_args()
    input_file = args.i
    output_file = args.o
    descr_file = args.s

    data = pd.read_csv(input_file, sep=CSV.get_delimiter(input_file), engine='python')
    sample_descr = pd.read_csv(descr_file, sep=CSV.get_delimiter(descr_file), engine='python')
    sample_name = sample_descr['Experiment'].iloc[0]

    data = data[columns]
    data = data.rename(columns={
        'Charge': 'Charges',
        'Intensity': sample_name + ' Intensity',
        'Peptide': 'Peptide Sequence',
        'Protein Start': 'Start',
        'Protein End': 'End'
    })
    data[sample_name + ' MaxLFQ Intensity'] = None
    data[sample_name + ' Match Type'] = 'MS/MS'

    data_spectral = data[['Peptide Sequence', 'Hyperscore']].groupby('Peptide Sequence').count().reset_index()
    data_spectral = data_spectral.rename(columns={'Hyperscore': sample_name + ' Spectral Count'})
    data = data.loc[data.groupby(['Peptide Sequence'])['Hyperscore'].idxmax()]
    data = data.drop(columns=['Hyperscore'])
    data = pd.merge(data, data_spectral, on=['Peptide Sequence'], how='left')

    data.to_csv(output_file, sep='\t', index=False)

    print('MSF_combined_peptide_maker.py: done')

    return None

# end of main()


if __name__ == '__main__':
    main()
