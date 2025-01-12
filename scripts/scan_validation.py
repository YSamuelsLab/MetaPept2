#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""
Created on Wed May 24 14:16:05 2023

@author: Bracha Erlanger Avigdor

"""

__version__ = '1.0.1'

import pathlib
import argparse
import os
import numpy as np
import pandas as pd


def is_valid_path(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return arg  # return arg as is


def get_scan_and_coverage_df(df, df_name):
    # extract scan number columns and pivot to make new dataframe
    # with [Sequence, Sample, and Scan number]
    scan_cols = ['Sequence'] + [x for x in df.columns if 'Scan number' in str(x)]
    scan_df = (df
               .query('Permutation_Index == 0')
               .filter(scan_cols)
               .melt(id_vars='Sequence',
                     value_vars=scan_cols,
                     var_name='Sample',
                     value_name='Scan_number').sort_values('Sequence'))

    # extract coverage columns and pivot to make new dataframe
    # with [Sequence, Sample, and coverage]
    target_cols = [x for x in df.columns if 'coverage' in str(x)]
    if not len(target_cols):  # the branch for MSFragger data
        target_cols = [x for x in df.columns if 'Hyperscore' in str(x)]
    if not len(target_cols):
        raise ValueError('ERROR: no columns to validate!')

    target_cols = ['Sequence'] + target_cols

    cov_df = (df
              .query('Permutation_Index == 0')
              .filter(target_cols)
              .melt(id_vars='Sequence',
                    value_vars=target_cols,
                    var_name='Sample_cov',
                    value_name='Coverage_' + df_name).sort_values('Sequence'))

    # merge the two dataframe to have [Sequence, Sample, Scan_number, Coverage]
    merge_df = (pd.merge(scan_df, cov_df, left_index=True, right_index=True)
                .drop(['Sequence_y', 'Sample_cov'], axis=1)
                .assign(Sample=lambda x: x['Sample'].str.replace('Scan number', ''))
                .rename(columns=dict(Sequence_x='Sequence'))
                .dropna(subset=['Scan_number', 'Coverage_' + df_name])
                ).reset_index(drop=True)

    return merge_df


def main(args):
    # read experiment IMP output
    df = pd.read_csv(args.experimental_file, engine='python')
    a_df = df.query("CDS.isna()", engine='python')

    # extract the cds only entries
    cds_df = df.query("CDS.notna() & nuORFs.isna()", engine='python')

    # remove sequences found in CDS (but keep nuORFs) database
    a_df = a_df.query("CDS.isna()", engine='python')

    A_df = get_scan_and_coverage_df(a_df, 'A')

    # read canonical  IMP output
    b_df = pd.read_csv(args.canonical_file, engine='python')
    B_df = get_scan_and_coverage_df(b_df, 'B')

    # merge experimental ad canonical on Sample and scan number to further compare the coverage
    A_B_merge_df = pd.merge(A_df, B_df, how='left', left_on=['Sample', 'Scan_number'],
                            right_on=['Sample', 'Scan_number'])

    # remove sequences that do not appear in the canonical database
    # select instances where:
    # -- scan numbers are the same
    # -- Sequences are different
    # -- coverage is higher or equal in the canonical vs experimental
    scan_number_list = (A_B_merge_df
                        .dropna(subset=['Sequence_y'])
                        .query('Sequence_x != Sequence_y and Coverage_A <= Coverage_B'))['Scan_number'].astype(
        np.int64).to_list()

    # search for scan numbers in the ABC dataframe and remove them
    # finally remove rows with no scan numbers
    cols = [x for x in a_df.columns if 'Scan number' in str(x)]
    a_df = (a_df
            .assign(**df.loc(axis=1)[cols]
                    .replace(scan_number_list, np.nan))
            .dropna(subset=cols, how='all')
            )

    # concatenate the validated experimental peptides
    # and the canonical (non-validated, original)
    scan_validation_df = pd.concat([a_df, cds_df])

    scan_validation_file = os.path.join(args.output_folder, 'IMP_scan_validation.csv')
    scan_validation_df.to_csv(scan_validation_file, index=False)

    print('scan validation file saved to: %s' % (scan_validation_file))


def make_parser():
    # define lambda function to validate file/path exists to use in add_argument
    is_valid = type = lambda x: is_valid_path(parser, x)

    parser = argparse.ArgumentParser(description="scan validation between experimental and canonical IMP outputs",
                                     epilog="Samuels Lab 2023", prog='')

    _input = parser.add_argument_group('input options')

    _input.add_argument('-a', '--experimental_file', metavar='', type=str, required=True,
                        help="experimental IMP filtered file")
    _input.add_argument('-b', '--canonical_file', metavar='', type=str, required=True,
                        help="canonical IMP filtered file")

    _output = parser.add_argument_group('output options')

    _output.add_argument('-o', '--output_folder', metavar='', type=is_valid, required=True,
                         help="folder location for output files (default: %(default)s)")

    parser.add_argument('-v', '--version', action='version', version="v%s" % (__version__))

    return parser


if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args()

    # call main
    main(args)
    print('Scan validation: done')
