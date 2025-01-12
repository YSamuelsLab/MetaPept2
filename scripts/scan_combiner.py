#!/usr/bin/env python3

import warnings
import sqlite3
import argparse
import pandas as pd
import re
import os
from itertools import islice
from CSVtools import CSV


def normalize_column_names(dataframe):
    warnings.simplefilter(action='ignore', category=FutureWarning)  # to suppress FutureWarning
    dataframe.columns = dataframe.columns.str.replace('[%()/*:]', '')
    dataframe.columns = dataframe.columns.str.strip().str.replace('[ .-]', '_')
    dataframe.columns = dataframe.columns.str.strip().str.replace('_+', '_')

    return dataframe

# end of normalize_column_names()


def seq_status(row, description):
    types = set()
    for name, val in row.items():
        if re.match(r'Scan_number_', name) and val and val > 0:
            sample, replica = islice(re.split(r'_+', name), 2, 4)
            desc = description[(description['Sample_Name'] == sample) & (description['Sample_Replica'].astype(str) == replica)]
            s_type = desc['Sample_Type'].iloc[0]
            types.add(s_type)

    types = ','.join(sorted(types))

    return types

# end of seq_status()


def combine(files, names, description_file, output, db_file=None):
    if len(files) != len(names):
        raise ValueError("The number of files must correspond to the number of experiment names")

    data_set = []
    col_names = ''
    delimiter = ','
    for i, file_name in enumerate(files):
        data = pd.read_csv(file_name, sep=CSV.get_delimiter(file_name), engine="python")
        normalize_column_names(data)

        if len(col_names):
            if col_names != ','.join(data.columns):
                raise ValueError("Combined files have different columns")
        else:
            col_names = ','.join(data.columns)

        data.insert(0, 'Exp', names[i])
        data_set.append(data)

    concat_data = pd.concat(data_set)

    if db_file:
        if os.path.exists(db_file):
            os.remove(db_file)
        connector = sqlite3.connect(db_file)
    else:
        connector = sqlite3.connect(':memory:')

    concat_data.to_sql(name='data', con=connector, index=False)
    cur = connector.cursor()

    sql = 'CREATE INDEX sequence_index ON data(Sequence);'
    cur.execute(sql)
    connector.commit()

    sql = 'CREATE TABLE ext_data AS SELECT t1.*, ' \
          '(SELECT GROUP_CONCAT(Exp) FROM (SELECT DISTINCT t2.Exp FROM data AS t2 ' \
          'WHERE t2.Sequence = t1.Sequence ORDER BY t2.Exp) AS t3) AS Experiment FROM data AS t1'
    cur.execute(sql)
    connector.commit()

    sql = 'CREATE INDEX ext_sequence_index ON ext_data(Sequence);'
    cur.execute(sql)
    connector.commit()

    db_created = False
    for name in names:
        if db_created:
            sql = 'INSERT OR IGNORE INTO combine SELECT * FROM ext_data WHERE Exp = "{}"'.format(name)
            cur.execute(sql)
        else:
            sql = 'CREATE TABLE combine AS SELECT * FROM ext_data WHERE Exp = "{}"'.format(name)
            cur.execute(sql)
            sql = 'CREATE UNIQUE INDEX combine_sequence_index ON combine(Sequence);'
            cur.execute(sql)
            db_created = True
        connector.commit()

    '''  # for future demands 
    sql = 'UPDATE combine AS t1 SET nuORFs = (SELECT GROUP_CONCAT(DISTINCT nuORFs) ' \
          'FROM ext_data AS t2 WHERE t2.Sequence = t1.Sequence)'
    cur.execute(sql)
    
    sql = 'UPDATE combine AS t1 SET CDS = (SELECT GROUP_CONCAT(DISTINCT CDS) ' \
          'FROM ext_data AS t2 WHERE t2.Sequence = t1.Sequence)'
    cur.execute(sql)
    '''

    sql = 'SELECT * FROM combine'
    combined_data = pd.read_sql(sql, connector)
    combined_data.drop(columns=['Exp'], inplace=True)

    desc_data = pd.read_csv(description_file, sep=CSV.get_delimiter(description_file), engine='python')
    combined_data['IMP_Status_over_sequence'] = combined_data.apply(seq_status, description=desc_data, axis=1)

    if True in combined_data.columns.str.contains('Hyperscore'):  # finding Best Hyperscores in MSFragger data
        columns = list(combined_data.columns) + ['BestHit_Hyperscore', 'BestHit_Deltascore', 'BestHit_MSFsample']
        combined_data = combined_data.reindex(columns=columns)

        def best_hyperscore(row):
            best_score = -1
            best_delta = None
            best_sample = None
            for col in row.index.to_list():
                if re.match('Hyperscore', col):
                    if row[col] > best_score:
                        best_score = row[col]
                        best_sample = re.sub('^Hyperscore[ _]*', '', col)
                        best_delta = row.filter(regex='^Delta[ _]+score[ _]*' + best_sample).iloc[0]

            row['BestHit_Hyperscore'] = best_score
            row['BestHit_Deltascore'] = best_delta
            row['BestHit_MSFsample'] = best_sample

            return row

        combined_data = combined_data.apply(best_hyperscore, axis=1)

    combined_data.to_csv(output, sep=delimiter, index=False)

    return combined_data

# end of combine()


def main():
    parser = argparse.ArgumentParser(description="IMP validation file combiner. It takes several IMP output files "
                                                 "and combine them together. The output file will have additional "
                                                 "column with experiment names.")

    parser.add_argument('-i', nargs='+', required=True, help="input file names")
    parser.add_argument('-n', nargs='+', required=True, help='list of experiment names (the order is important)')
    parser.add_argument('-s', default='sample_description.csv', required=True,
                        help='tab delimited file with the description of samples: '
                             'Source_File	Sample_Name	Sample_Replica	Sample_Type Group')
    parser.add_argument('-o', default='combined_file.csv', required=False, help='output file with joined data')
    parser.add_argument('-db', required=False, help='database file (if not specified, the data will be stored in memory')

    args = parser.parse_args()

    input_file_list = args.i
    experiment_list = args.n
    description_file = args.s
    output_file = args.o
    db_file = args.db

    try:
        combine(input_file_list, experiment_list, description_file, output_file, db_file)
    except Exception as err:
        print("Something went wrong: {}".format(err))

    print('Scan combiner: done')

# end of main()


if __name__ == '__main__':
    main()
