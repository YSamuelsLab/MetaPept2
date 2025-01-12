#!/usr/bin/env python3

import re
import os
import argparse
import pandas as pd
import sqlite3
import warnings
import glob
import csv
import gzip

# INFO: some peptides in PEAKS output files have no Intensity value!
# In that case you can find a non-zero replica count but zero for all Replica Intensities and zero for Intensity Sum.

# You need a sample description file to run the script:

# ######### Sample description file format (tab delimited) #####################################
# Source_File                     Sample_Name  Sample_Replica  Sample_Type  Group Experiment   #
# HF2_YS_14535_1043_18082021.raw  SCC1         R1              ko           1     SCC1 1 HLA-I #
# HF2_YS_14535_1042_18082021.raw  SCC1         R2              ko           1     SCC1 1 HLA-I #
# HF2_YS_14535_1047_18082021.raw  WT           R1              wt           1     SCC1 1 HLA-I #
# ##############################################################################################

# Replicas with the same Group name will be processed together and saved in the individual output file

# Annotated files must contain category names in their paths or names (the position is not important)
# For example: frameshift.sample_name.csv.gz.pep.annotated.csv.gz or sample_name.prio2.csv.gz.pep.annotated.csv.gz
# The default categories are: frameshift, prio1, prio2, prio3

# the pattern for human and mouse HLA allele variants (type I and type II)
HLA_pattern = r'HLA[-_][ABCEG]\d|D[PQR][AB]\d+[-_]|H[-_]?2[-_][DKLQI]'


def get_delimiter(file_path):
    sniffer = csv.Sniffer()
    if re.search(r'.gz$', file_path):
        with gzip.open(file_path, "rt") as f:
            header = f.readline()
    else:
        with open(file_path, "rt") as f:
            header = f.readline()

    delimiter = sniffer.sniff(header).delimiter

    if delimiter not in ['\t', ',']:
        file_name = re.sub(r'.*/', '', file_path)
        raise ValueError('Can not determine CSV delimiter for {} file'.format(file_name))

    return delimiter

# end of get_delimiter()


def counter():  # just for debug
    count = 0

    def item():
        nonlocal count
        count += 1
        return count

    return item

# end of counter()


def get_category(file_path, cat_aliases):
    category = None
    for cat in cat_aliases:
        if file_path.find(cat) > -1:
            if category:
                raise ValueError("ERROR: incorrect category in the file name or path {}".format(file_path))
            else:
                category = cat
    if not category:
        raise ValueError("ERROR: no category in the file name or path {}".format(file_path))

    return category

# end of get_category()


def drop_columns(dataframe, column_list):
    for col in column_list:
        if col in dataframe.columns:
            dataframe.drop(col, axis=1, inplace=True)

    return dataframe

# end of drop_columns()


def normalize_column_names(dataframe):
    warnings.simplefilter(action='ignore', category=FutureWarning)  # to suppress FutureWarning
    dataframe.columns = dataframe.columns.str.replace('[%()/*:]', '')
    dataframe.columns = dataframe.columns.str.strip().str.replace('[ .-]', '_')
    dataframe.columns = dataframe.columns.str.strip().str.replace('_+', '_')

    return dataframe

# end of normalize_column_names()


def normalize_allele_names(dataframe):
    warnings.simplefilter(action='ignore', category=FutureWarning)  # to suppress FutureWarning
    dataframe['HLA_allele'] = dataframe['HLA_allele'].str.replace('^.$', '')
    dataframe['HLA_allele'] = dataframe['HLA_allele'].str.replace('[*:]', '')
    dataframe['HLA_allele'] = dataframe['HLA_allele'].str.replace('-', '_')

    return dataframe

# end of normalize_allele_names()


def transform_allele_columns(dataframe):
    hla_alleles = []
    for col in dataframe.columns:
        if re.search(HLA_pattern, col):
            hla_alleles.append(col)

    for i, allele in enumerate(hla_alleles):
        dataframe.rename(columns={allele: 'netMHC_rank_' + str(i+1)}, inplace=True)
        dataframe.insert(len(dataframe.columns), 'HLA_allele_' + str(i+1), allele)

    return hla_alleles

# end of transform_allele_columns()


def check_allele_consistency(dataframe):
    allele_names = []
    for col in dataframe.columns:
        if re.search(r'HLA_allele_.*\d', col):
            hla_alleles = list(dataframe[col].unique())
            if len(hla_alleles) == 1:
                allele_names.append(hla_alleles[0])
            else:
                return []

    return allele_names

# end of check_allele_consistency()


def replace_hla_column_names(column_names, dataframe):
    i_column = 0
    for col in dataframe.columns:
        if re.search(r'netMHC_rank_.*\d', col):
            if i_column < len(column_names):
                dataframe.rename(columns={col: column_names[i_column]}, inplace=True)
            else:
                return False
            i_column += 1

    for col in dataframe.columns:
        if re.search(r'HLA_allele_.*\d', col):
            dataframe.drop(col, axis=1, inplace=True)

    return True

# end of replace_hla_column_names()


def combine(input_dir, sample_file, db_file, q_threshold, alc_threshold, rank_threshold, output_file, decoy, cat_aliases):
    print('Data preparing ...', end='', flush=True)
    all_data = []
    for file_path in glob.glob(input_dir + '/**/*.pep.annotated.csv.gz', recursive=True):
        category = get_category(file_path, cat_aliases)
        data = pd.read_csv(file_path, compression='gzip', sep=get_delimiter(file_path), engine='python')

        normalize_column_names(data)
        normalize_allele_names(data)
        drop_columns(data, ['Denovo_score', 'Predict_RT'])  # drop some columns because of the different file formats
        data = data[data['Decoy'] == 'D'] if decoy else data[data['Decoy'] != 'D']
        data.insert(len(data.columns), 'Databases_PRISM', category)
        transform_allele_columns(data)  # it needs to unify table for samples with different HLA alleles
        all_data.append(data)

    if not len(all_data):
        raise ValueError('ERROR: no valid files in the input directory!')
    all_data = pd.concat(all_data)

    if db_file:
        if os.path.exists(db_file):
            os.remove(db_file)
        connector = sqlite3.connect(db_file)
    else:
        connector = sqlite3.connect(':memory:')

    sample_description = pd.read_csv(sample_file, sep=get_delimiter(sample_file), quotechar='"', low_memory=False)
    sample_description.to_sql(name='description', con=connector, index=False)

    all_data.index += 1  # to start the index with 1, but not zero
    all_data.index.name = 'Rec_ID'

#   print(list(all_data.columns))  # for debug
    all_data.to_sql(name='data', con=connector)

    cur = connector.cursor()
    sql = 'CREATE INDEX source_file_index ON description(Source_File);'
    cur.execute(sql)

    sql = 'CREATE TABLE ext_data AS SELECT * FROM data INNER JOIN description USING(Source_File);'
    cur.execute(sql)
    connector.commit()

    sql = 'CREATE INDEX group_index ON ext_data("Group");'
    cur.execute(sql)

    sql = 'SELECT COUNT(*) AS n FROM ext_data;'
    n_rec_ext_data = pd.read_sql(sql, connector).iloc[0]['n']
    sql = 'SELECT COUNT(*) AS n FROM data;'
    n_rec_data = pd.read_sql(sql, connector).iloc[0]['n']

    print('OK')
    print('Combining groups:')

    if n_rec_ext_data != n_rec_data:
        print("You probably have an incorrect file with the description of the samples")
        sql = 'SELECT DISTINCT Source_File AS Files FROM data;'
        data_files = '\n'.join(pd.read_sql(sql, connector)['Files'].tolist())
        print('\nfiles in the data:')
        print(data_files)

        sql = 'SELECT DISTINCT Source_File AS Files FROM description;'
        desc_files = pd.read_sql(sql, connector)['Files'].tolist()
        print('\nfiles in the description:')
        print(desc_files)

        exit(1)

    # FROM Andreas's emails:
    # I typically filter by Q (e.g. Q<0.01 = 1%FDR) and by NetMHC predictions (e.g. rank<2%).
    # I generate this column in Spotfire. Best Q is the lowest Q value for a peptide sequence among all merged samples.
    # I.e., when a peptide is reliably identified in sample 1 (e.g. Q<0.01) and in sample 2 Q>0.01,
    # the peptide is still counted as identified in both samples. This increases the overlap between merged samples.
    #
    # "best ALC" gives you the highest ALC of a sequence over all combined analyses.
    # If you compare "best ALC" and "best Q" you will see that the values correlate.
    # Typically, the sequence with the best ALC will have the lowest best Q value. However, there can be exceptions.

    sql = 'SELECT DISTINCT "Group" FROM ext_data;'  # DON'T REMOVE DOUBLE QUOTES AROUND "Group"!!!
    groups = pd.read_sql(sql, connector)
    group_number = len(groups['Group'])
    for group in groups['Group']:
        sql = 'DROP TABLE IF EXISTS tmp_data;'
        cur.execute(sql)
        connector.commit()

        sql = 'SELECT * FROM ext_data WHERE "Group" = "{}";'.format(group)
        group_data = pd.read_sql(sql, connector)
        hla_allele_names = check_allele_consistency(group_data)
        if not len(hla_allele_names):
            raise ValueError('ERROR: the group {} has inconsistent set of alleles'.format(group))

        if not replace_hla_column_names(hla_allele_names, group_data):
            raise ValueError('ERROR: column names were not replaced for the group {}'.format(group))

        group_data.to_sql(name='tmp_data', con=connector)
        sql = 'CREATE INDEX sequence_index ON tmp_data(Sequence);'
        cur.execute(sql)
        connector.commit()

        sql = 'SELECT DISTINCT Sequence FROM tmp_data;'
        peptides = pd.read_sql(sql, connector)
        data_output = []
        # n_debug = counter()  # for debug
        for peptide in peptides['Sequence']:
            # if n_debug() > 10:  # for debug
            #     break

            print('group {}  {}'.format(group, peptide) + ' ' * 40, end="\r")
            sql = 'SELECT * FROM tmp_data WHERE Sequence = "{}" ORDER BY Q, netMHC_rank LIMIT 1;'.format(peptide)
            best_peptide_rec = pd.read_sql(sql, connector)

            if len(best_peptide_rec.index) and best_peptide_rec.iloc[0]['Q'] < q_threshold:
                best_q = best_peptide_rec.iloc[0]['Q']  # minimum false discovery rate (FDR)
                best_sample = best_peptide_rec.iloc[0]['Sample_Name']
                best_replica = best_peptide_rec.iloc[0]['Sample_Replica']

                sql = 'SELECT ALC FROM tmp_data WHERE Sequence = "{}" ORDER BY ALC DESC LIMIT 1;'.format(peptide)
                best_alc_rec = pd.read_sql(sql, connector)
                best_alc = best_alc_rec.iloc[0]['ALC']  # average local confidence (ALC)

                # this is the legacy code for rank filtering, the data cleanup step was added below
                best_hla = best_peptide_rec.iloc[0]['HLA_allele']
                filtered_hla = best_hla if re.search(HLA_pattern, best_hla) and \
                                           best_hla in best_peptide_rec.columns and \
                                           best_peptide_rec.iloc[0][best_hla] < rank_threshold else ''

                sql = 'SELECT group_concat(DISTINCT Category) AS Categories FROM tmp_data WHERE Sequence = "{}";'.format(peptide)
                categories_rec = pd.read_sql(sql, connector)
                categories = categories_rec.iloc[0]['Categories']

                sql = 'SELECT group_concat(DISTINCT Sample_Type) AS Types FROM tmp_data ' \
                      'WHERE Sequence = "{}";'.format(peptide)
                status_rec = pd.read_sql(sql, connector)
                status = status_rec.iloc[0]['Types']

                sql = 'SELECT group_concat(DISTINCT Location) AS Location FROM tmp_data ' \
                      'WHERE Sequence = "{}";'.format(peptide)
                status_rec = pd.read_sql(sql, connector)
                location = status_rec.iloc[0]['Location']

                sql = 'SELECT group_concat(DISTINCT Gene) AS Gene FROM tmp_data ' \
                      'WHERE Sequence = "{}";'.format(peptide)
                status_rec = pd.read_sql(sql, connector)
                gene = status_rec.iloc[0]['Gene']

                sql = 'SELECT group_concat(DISTINCT Symbol) AS Symbol FROM tmp_data ' \
                      'WHERE Sequence = "{}";'.format(peptide)
                status_rec = pd.read_sql(sql, connector)
                symbol = status_rec.iloc[0]['Symbol']

                sql = 'SELECT group_concat(DISTINCT ORF_location) AS ORF_location FROM tmp_data ' \
                      'WHERE Sequence = "{}";'.format(peptide)
                status_rec = pd.read_sql(sql, connector)
                orf_location = status_rec.iloc[0]['ORF_location']

                sql = 'SELECT group_concat(Databases_PRISM) AS Databases_PRISMs FROM ' \
                      '(SELECT DISTINCT Databases_PRISM FROM tmp_data WHERE Sequence = "{}" ORDER BY Databases_PRISM);'.format(peptide)
                db_prism_rec = pd.read_sql(sql, connector)
                db_prism = db_prism_rec.iloc[0]['Databases_PRISMs']

                sql = 'SELECT group_concat(Sample_Name) AS Sample_Names FROM ' \
                      '(SELECT DISTINCT Sample_Name FROM tmp_data WHERE Sequence = "{}" ORDER BY Sample_Name);'.format(peptide)
                sample_rec = pd.read_sql(sql, connector)
                samples = sample_rec.iloc[0]['Sample_Names']

                sql = 'SELECT group_concat(scan) AS Scans FROM ' \
                      '(SELECT DISTINCT Sample_Name || "/" || Sample_Replica || "=" || Scan AS scan ' \
                      'FROM tmp_data WHERE Sequence = "{}" ' \
                      'ORDER BY Sample_Name, Sample_Replica);'.format(peptide)
                scans_rec = pd.read_sql(sql, connector)
                all_scans = scans_rec.iloc[0]['Scans']

                sql = 'SELECT DISTINCT Sample_Name FROM description WHERE "Group" = "{}" ' \
                      'ORDER BY Sample_Name;'.format(group)
                samples_rec = pd.read_sql(sql, connector)
                intensity_sum = 0
                intensities = []
                replica_counts = []
                for sam_index, sam_row in samples_rec.iterrows():
                    sample_name = sam_row['Sample_Name']
                    sql = 'SELECT DISTINCT Sample_Replica FROM tmp_data ' \
                          'WHERE Sequence = "{}" AND Sample_Name = "{}";'.format(peptide, sample_name)
                    replica_count_rec = pd.read_sql(sql, connector)
                    replica_counts.append({'sample': sample_name,
                                           'rep_count': len(replica_count_rec['Sample_Replica'].values)
                                           })

                    sql = 'SELECT DISTINCT Sample_Replica FROM description ' \
                          'WHERE Sample_Name = "{}" ORDER BY Sample_Replica;'.format(sample_name)
                    replica_rec = pd.read_sql(sql, connector)
                    for rep_index, rep_row in replica_rec.iterrows():
                        sample_replica = rep_row['Sample_Replica']
                        sql = 'SELECT Intensity FROM tmp_data ' \
                              'WHERE Sequence = "{}" AND Sample_Name = "{}" AND Sample_Replica = "{}" ' \
                              'ORDER BY Intensity DESC LIMIT 1;'.format(peptide, sample_name, sample_replica)
                        intensity_rec = pd.read_sql(sql, connector)
                        intensity = intensity_rec.iloc[0]['Intensity'] \
                            if len(intensity_rec.index) and intensity_rec.iloc[0]['Intensity'] else 0

                        intensities.append({'replica': '_'.join([str(sample_name), str(sample_replica)]),
                                            'intensity': intensity})
                        intensity_sum += intensity

                hla_columns = [col for col in best_peptide_rec.columns if re.match(HLA_pattern, col)]
                base_columns = ['Source_File', 'Feature', 'Scan', 'ALC', 'Length', 'RT', 'Mass', 'ppm', 'ID',
                                'Location_count', 'Genome', 'Location', 'Sequence', 'Top_location_count',
                                'Top_location_count_no_decoy', 'Q', 'Gene', 'Symbol', 'ORF_location', 'HLA_allele',
                                'netMHC_rank']
                data_output_rec = best_peptide_rec[base_columns + hla_columns]
                data_output_rec = data_output_rec.assign(
                    Filtered_HLA_allele=filtered_hla,
                    Best_Q=best_q,
                    Best_ALC=best_alc,
                    Best_Q_replica='/'.join([str(best_sample), str(best_replica)]),
                    Categories=categories,
                    Status_over_sequence=status,
                    Databases_PRISM=db_prism,
                    Samples=samples,
                    All_scans=all_scans
                )

                # Updates for some fields to integrate all found entries
                data_output_rec.loc[0, 'Location'] = location
                data_output_rec.loc[0, 'Gene'] = gene
                data_output_rec.loc[0, 'Symbol'] = symbol
                data_output_rec.loc[0, 'ORF_location'] = orf_location

                replica_report = pd.DataFrame(
                    [[item['rep_count'] for item in replica_counts]],
                    columns=[item['sample'] for item in replica_counts]
                )
                normalize_column_names(replica_report)
                data_output_rec = pd.concat([data_output_rec, replica_report], axis=1)

                data_output_rec.insert(len(data_output_rec.columns), 'Intensity_Sum', intensity_sum)

                intensity_report = pd.DataFrame(
                    [[item['intensity'] for item in intensities]],
                    columns=['Intensity_' + item['replica'] for item in intensities]
                )
                normalize_column_names(intensity_report)
                data_output_rec = pd.concat([data_output_rec, intensity_report], axis=1)
                data_output.append(data_output_rec)

        data_output = pd.concat(data_output)
        data_output = data_output.replace(['-'], '')
        path = re.split(r'/', output_file)
        if group_number > 1:  # if there is only one group the output file name won't be modified
            path[-1] = group + '_' + path[-1]
        group_output_file = '/'.join(path)

        data_output = data_output[data_output['Best_ALC'] >= alc_threshold]
        data_output = data_output[data_output['netMHC_rank'] < rank_threshold]

        if re.search(r'\.gz$', group_output_file):
            data_output.to_csv(group_output_file, compression='gzip', sep='\t', index=False)
        else:
            data_output.to_csv(group_output_file, sep='\t', index=False)

        print('group {} ...Ok'.format(group) + ' ' * 50)

# end of combine()


def main():
    parser = argparse.ArgumentParser(description="PRISM combiner. It takes all {category}.*.pep.annotated.csv.gz files "
                                                 "and combines them into one output file according to "
                                                 "the Sample Description file.")

    parser.add_argument('-i', required=True, help="input directory (all combined files must be in the directory)")
    parser.add_argument('-s', default='sample_description.csv', required=True,
                        help='tab delimited file with the description of samples: '
                             'Source_File	Sample_Name	Sample_Replica	Sample_Type Group')
    parser.add_argument('-db', required=False,
                        help='database file (if not specified, the data will be stored in memory)')
    parser.add_argument('-q', default=100, type=float, required=False, help='Q threshold: 0.1 = 10%%FDR (default: no filtering)')  # don't remove double "%"
    parser.add_argument('-best_alc', default=-100, type=float, required=False, help='Best ALC threshold (default: no filtering)')
    parser.add_argument('-r', default=2.0, type=float, required=False, help='Rank threshold: SB < 0.5; WB < 2.0 (default: 2.0)')
    parser.add_argument('-o', default='combined_results.csv.gz', required=False, help='output file with results')
    parser.add_argument('-cat', default=['frameshift', 'prio1', 'prio2', 'prio3'], nargs='+', required=False, help='category aliases for running PRISM')
    parser.add_argument('-D', action='store_true', help='combine decoy peptides')

    args = parser.parse_args()

    input_dir = re.sub(r'/$', '', args.i)
    sample_file = args.s
    db_file = args.db
    q_threshold = args.q
    alc_threshold = args.best_alc
    rank_threshold = args.r
    output_file = args.o
    decoy = args.D
    cat_aliases = args.cat

    combine(input_dir, sample_file, db_file, q_threshold, alc_threshold, rank_threshold, output_file, decoy, cat_aliases)

    print("PRISM combiner: done")

# end of main()


if __name__ == '__main__':
    main()
