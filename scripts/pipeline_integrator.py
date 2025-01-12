#!/usr/bin/env python3

import argparse
import re
import pandas as pd
import warnings
import Levenshtein
from CSVtools import CSV

MQ_COLUMNS2REMOVE = [r'^Charge_.*', r'^Mass_.*', r'term_cleavage_window']
PRISM_COLUMNS2REMOVE = [r'^Location_count$', r'^Genome$', r'^Top_location_count_no_decoy$',
                        r'^netMHC_rank$', r'^HLA_.*', r'^Intensity_.*', r'^D[PQR][AB]\d+_', r'^H[-_]?2[-_][DKLQI]']


def normalize_column_names(dataframe):
    warnings.simplefilter(action='ignore', category=FutureWarning)  # to suppress FutureWarning
    dataframe.columns = dataframe.columns.str.replace('[%()/*:]', '')
    dataframe.columns = dataframe.columns.str.strip().str.replace('[ .-]', '_')
    dataframe.columns = dataframe.columns.str.strip().str.replace('_+', '_')

    return dataframe


def delete_columns(patterns, data):
    n = 0
    for pattern in patterns:
        if isinstance(data, pd.DataFrame):
            columns = [col for col in data.columns if re.search(pattern, col)]
            if len(columns):
                data.drop(columns=columns, inplace=True)
                n += len(columns)
        elif isinstance(data, pd.Series):
            columns = [col for col in data.index if re.match(pattern, col)]
            if len(columns):
                data.drop(index=columns, inplace=True)
                n += len(columns)

    return n


def fuzzy_map_scan_column_names(raw2replica, mq):
    scan_columns = list(mq.filter(regex='Scan.*number').columns)

    if len(raw2replica.keys()) != len(scan_columns):
        raise ValueError('Table comparison error: '
                         '{} raw files in Sample description file '
                         'but only {} `Scan number` columns.'.format(len(raw2replica.keys()), len(scan_columns)))

    replica2column = {}
    for replica in raw2replica.keys():
        replica_score = 0
        for col_name in scan_columns:
            similarity_score = Levenshtein.ratio(replica, col_name)
            if replica_score < similarity_score:
                replica_score = similarity_score
                replica2column[replica] = col_name

    return replica2column


def map_scan_column_names(replica2experiment, mq):  # this ugly code is a legacy of the pipeline transformation
    scan_columns = list(mq.filter(regex='^Scan_number').columns)

    if len(replica2experiment.keys()) != len(scan_columns):
        raise ValueError('Table comparison error: '
                         '{} raw files in Sample description file '
                         'but only {} `Scan number` columns.'.format(len(replica2experiment.keys()), len(scan_columns)))

    replica2column = {}
    for replica, experiment in replica2experiment.items():
        column_name = pd.DataFrame([], columns=['Scan number ' + experiment.strip()])
        norm_column_name = normalize_column_names(column_name).columns[0]  # create the same column format

        column_not_found = True
        for col_name in scan_columns:
            if norm_column_name == col_name:
                replica2column[replica] = col_name
                column_not_found = False
                break
        if column_not_found:
            raise ValueError('The replica {} was not found in IMP file! \n'
                             'The file columns: {}'.format(replica, ','.join(scan_columns)))

    return replica2column


def get_filenames(prefix_set, name):
    new_names = []
    for prefix in prefix_set:
        path = re.split(r'/', name)
        path[-1] = prefix + '_' + path[-1]
        new_names.append('/'.join(path))

    return new_names


def combine(description_file, mq_file, prism_file, output_file, strict_mode=False):
    try:
        mq = pd.read_csv(mq_file, sep=CSV.get_delimiter(mq_file), engine='python')
        prism = pd.read_csv(prism_file, sep=CSV.get_delimiter(prism_file), engine='python')
        description = pd.read_csv(description_file, sep=CSV.get_delimiter(description_file), engine='python')
    except Exception as err:
        print(err)
    else:
        normalize_column_names(mq)
        normalize_column_names(prism)
        mq4unique = mq.copy()  # to keep the original data structure for MQ-unique
        delete_columns(MQ_COLUMNS2REMOVE, mq)
        description['Sample_Replica'] = description['Sample_Replica'].astype(str)
        description['Sample_Name'] = description['Sample_Name'].astype(str)

        replica2raw = dict(
            zip(description['Sample_Name'] + '_' + description['Sample_Replica'], description['Source_File'])
        )
        # replica2column_name_fuzzy = fuzzy_map_scan_column_names(replica2raw, mq)

        # the fuzzy match is replaced with a strict match
        replica2experiment = dict(
            zip(description['Sample_Name'] + '_' + description['Sample_Replica'], description['Experiment'])
        )
        replica2column_name = map_scan_column_names(replica2experiment, mq)

        raw2replica = dict([(value, key) for key, value in replica2raw.items()])

        # to delete columns from IMP table with the same name in PRISM table
        column_intersection = set(prism.columns).intersection(mq.columns)

        shared_rows = []
        prism_unique_rows = []
        # debug_n = 0
        for prism_index, prism_row in prism.iterrows():
            prism_peptide = prism_row['Sequence']
            prism_scan_number = prism_row['Scan']
            if prism_row['Source_File'] in raw2replica and raw2replica[prism_row['Source_File']] in replica2column_name:
                replica = raw2replica[prism_row['Source_File']]
                md_scan_colum = replica2column_name[replica]
                if strict_mode:
                    # IMP and DENOVO have the same scan number for the peptide in the replica
                    mq_pept_rows = mq.loc[(mq['Sequence'] == prism_peptide) & (mq[md_scan_colum] == prism_scan_number)]
                else:
                    # IMP and DENOVO can have different scan numbers for the peptide in the replica
                    mq_pept_rows = mq.loc[(mq['Sequence'] == prism_peptide) & (mq[md_scan_colum] > 0)]
                if len(mq_pept_rows.index):
                    delete_columns(PRISM_COLUMNS2REMOVE, prism_row)
                    for mq_index, mq_pept_row in mq_pept_rows.iterrows():
                        concat = pd.concat([prism_row, mq_pept_row.drop(labels=column_intersection)])
                        shared_rows.append(pd.Series([replica], index=['Best_PRISM_Replica']).append(concat))

                else:
                    prism_unique_rows.append(prism_row)
            else:
                raise ValueError("Can not find DENOVO's sample name in IMP table")

        shared_data_output = pd.DataFrame(shared_rows)
        prism_unique_data_output = pd.DataFrame(prism_unique_rows)

        shared_sequences = list(shared_data_output['Sequence'])
        mq_unique_data_output = mq4unique[~mq4unique['Sequence'].isin(shared_sequences)].copy()
        delete_columns(MQ_COLUMNS2REMOVE, mq_unique_data_output)

        combined_output_file, prism_unique_output_file, mq_unique_output_file = get_filenames(
            ['combined', 'denovo_unique', 'imp_unique'], output_file)
        if re.search(r'\.gz$', output_file):
            shared_data_output.to_csv(combined_output_file, compression='gzip', sep='\t', index=False)
            prism_unique_data_output.to_csv(prism_unique_output_file, compression='gzip', sep='\t', index=False)
            mq_unique_data_output.to_csv(mq_unique_output_file, compression='gzip', sep='\t', index=False)
        else:
            shared_data_output.to_csv(combined_output_file, sep='\t', index=False)
            prism_unique_data_output.to_csv(prism_unique_output_file, sep='\t', index=False)
            mq_unique_data_output.to_csv(mq_unique_output_file, sep='\t', index=False)

    return None


# end of combine()


def main():
    parser = argparse.ArgumentParser(description='A script for the scan integration of DENOVO and IMP pipelines')
    parser.add_argument('-imp', required=True, help='IMP file (MaxQuant/MSFragger)')
    parser.add_argument('-denovo', required=True, help='DENOVO file (PRISM)')
    parser.add_argument('-s', required=True, help='Sample Description file')
    parser.add_argument('-o', default='scan_integration.csv', required=False, help='output file')
    parser.add_argument('-strict', action='store_true', help='strict mode: peptide scan numbers have to be matched '
                                                             'between IMP and DENOVO outputs')

    args = parser.parse_args()
    imp_file = args.imp
    denovo_file = args.denovo
    description_file = args.s
    output_file = args.o
    strict = args.strict

    try:
        combine(description_file, imp_file, denovo_file, output_file, strict)
    except Exception as err:
        print("Something went wrong:", err)

    print('Integration: done')


# end of main()


if __name__ == '__main__':
    main()
