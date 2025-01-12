#!/usr/bin/env python3

import argparse
import re
import csv
import sys
import gzip
import pandas as pd
import numpy as np
from difflib import SequenceMatcher
from CSVtools import CSV


COLUMNS_SET1 = ['^Seq$', '^Chimera$', '^Out_Frame$', '^Pattern$', '^DB$', 'Source_db$', '^Source_File$',
           '^Feature$', '^Scan$', '^Length$', '^RT$', '^Location$', '^Sequence$', '^Top_location_count$',
           '^Gene$', '^Symbol$', '^ORF_location$', '^Best_Q$', '^Best_ALC$', '^Categories$', '^Status_over_sequence$',
           '^IMP_Status_over_sequence$', '^Databases_PRISM$', '^Samples$']

# Here is the place for sample specific columns

COLUMNS_SET2 = ['^Gene_ID$', '^Gene_Name$', '^nuORFs$', '^Hits$', '^Proteins$', '^Leading_razor_protein$', '^Score$', '^Experiment_',
           '^Intensity_', '^HLA_[^_]+$', '^HLA_[^_]+_AffnM', '^H_2_[^_]+$', '^H_2_[^_]+_AffnM', '^Scan_number_', '^Experiment$', '^Best_coverage$', '^Best_delta_score$',
           '^BestHit_Hyperscore$', '^BestHit_Deltascore$', '^Integration$',
           '^HLA_atlas_WEB$', '^HLA_atlas_PRISM$', '^IEDB$', '^iREAD$', '^GATK$', '^Extra$']


def get_samples(data):
    samples = []
    for column in list(data.columns):
        search_exp = re.search(r'Scan_number_(.*)_(?:HLA|MHC)_.*', column)
        if search_exp:
            sample, replica = re.split('_', search_exp.group(1))
            if sample not in samples:
                samples.append(sample)

    return samples

# end of get_samples()


def walk_across_samples(row, sample, update_column):
    num_replicas = 0
    for column in list(row.index):
        if re.search(r'Scan_number_{}_.*_(?:HLA|MHC)_.*'.format(sample), column):
            if not np.isnan(row[column]):
                num_replicas += 1

    row[update_column] = num_replicas
    return row

# end of walk_across_samples()


def parse_header(body):  # to parse Extra column
    records = {}
    if not pd.isnull(body):
        num = 0
        for item in re.split(',', body):
            peptide, header = re.split(':', item)
            peptide, pos = re.split('_', peptide)
            pos, pseq_len = re.split('/', pos)

            # FORMAT (DEFAULT): AAAAGALPR_8/17:P1_emptyA_count3_F487_ENST00000611612.2_TTC34_LPRP_GGCCTTTGCC
            # FORMAT (OOF): QATTVLHIL_132/180:P1_Asite_count1_F52_ENST00000629496.3_DDX3X_PPHLRNREATKGF_TIKTVQGGVLAKIR_AGGTTTCTAC
            site_type, rest = re.split('_count', header)
            items = re.split('_', rest)
            if len(items) == 6:  # FORMAT DEFAULT
                count, f_num, gene_id, gene_name, out_frame, pattern = items
                in_frame = None
            elif len(items) == 7:  # FORMAT OOF
                count, f_num, gene_id, gene_name, in_frame, out_frame, pattern = items
            else:
                raise ValueError('ERROR: wrong FASTA header')

            if peptide not in records:
                records[peptide] = []

            source_db = None
            if re.search('emptyA', site_type):
                source_db = 'Aeffect'
            elif re.search('Psite|Asite', site_type):
                source_db = 'Peffect'

            records[peptide].append({
                'num': num + 1,
                'pos': pos,
                'pseq_len': pseq_len,
                'site_type': site_type,
                'db': site_type,
                'source_db': source_db,
                'pattern': pattern,
                'in_frame': in_frame,
                'out_frame': out_frame,
                'gene_id': gene_id,
                'gene_name': gene_name
            })

    return records

# end of parse_header()


def matching(string1, string2):
    # This function is to write OOF part of peptide in lower character
    # The function was taken from 'pip_peptidomics_analysis_aberrant_Deborah.py'
    # The function was replaced with overlap() function due to bugs in PRISM (some canonical peptides reported as non-canonical)

    if string1 in string2:
        match = string1.lower()
    else:

        match = SequenceMatcher(None, string1[::-1], string2[::-1]).find_longest_match(0, len(string1), 0, len(string2))
        match = string1[::-1][0: match.size].lower() + string1[::-1][match.size:].upper()
        match = match[::-1]

    return match

# end of matching()


def overlap(peptide, pos, pseq_len, off_frame, in_frame=None):
    # This function was written to fix problems with `pip_peptidomics_analysis_aberrant_Deborah.py' script

    if in_frame is not None:  # OOF database
        break_pos = len(in_frame) - pos + 1
        if break_pos < 0:
            break_pos = 0
        pep_head = peptide[:break_pos].upper()
        pep_tail = peptide[break_pos:].lower()
    else:  # Default
        break_pos = pseq_len - len(off_frame) - pos + 1
        if break_pos < 0:
            break_pos = 0
            if peptide not in off_frame:
                # we do not check it for OOF because OFFFrame sequence in FASTA file can be much longer than
                # `off_frame` value. But for the default database we do it.
                raise ValueError('Can not find the peptide {} in the off frame sequence {}'.format(peptide, off_frame))

        pep_head = peptide[:break_pos].upper()
        pep_tail = peptide[break_pos:].lower()
        if (len(pep_head) and not re.match(pep_tail, off_frame.lower())) or \
                (not len(pep_head) and pep_tail not in off_frame.lower()):
            raise ValueError('Wrong peptide tail {} for the off frame sequence {}'.format(pep_tail, off_frame))

    chimera = pep_head + pep_tail

    if len(chimera) != len(peptide):
        raise ValueError('Wrong chimera sequence {} ({}/{})'.format(chimera, peptide, off_frame))

    return chimera

# end of matching()


def multi_matching(pseq, extra, unique=False):
    # Wrap-up function for overlap()

    matches = []
    if pd.notna(pseq) and pd.notna(extra):
        header_data = parse_header(extra)
        if pseq not in header_data:
            raise ValueError('ERROR: can not find the peptide sequence ({}) in the FASTA file header'.format(pseq))

        for peptide in header_data:
            for rec in header_data[peptide]:
                matches.append(overlap(peptide, int(rec['pos']), int(rec['pseq_len']), rec['out_frame'], rec['in_frame']))

        if unique:
            matches = sorted(set(matches))

    return ','.join(matches) if len(matches) else None

# end of multi_matching()


def set_from_extra(pseq, body, category, field, unique=False):
    # to create new columns based on Extra column information

    records = parse_header(body)
    patterns = []
    if re.search('Extra', str(category)):  # only `Extra` category will be treated
        if pseq in records:
            for item in records[pseq]:
                patterns.append(item[field])
        else:
            return None

        if unique:
            patterns = sorted(set(patterns), key=str.casefold)

        return ','.join(patterns) if len(patterns) else None

    return None

# end of set_from_extra()


def make_table(file_name, sb_threshold, wb_threshold, all_data=False):
    csv.field_size_limit(sys.maxsize)
    data = pd.read_csv(file_name, sep=CSV.get_delimiter(file_name), engine='python')
    # remove canonical to mitigate Denovo (PRISM) annotation errors
    data = data[data['Canonical'].isna()]

    # column set
    column_names = COLUMNS_SET1

    # combine MQ sample information for each peptide: how many replicas have exposed the peptide
    mq_column_names = []
    for sample in get_samples(data):
        mq_column_name = 'IMP_' + sample
        data = data.assign(mq_column_name=0)
        mq_column_names.append('^' + mq_column_name + '$')
        data = data.apply(walk_across_samples, args=(sample, mq_column_name), axis=1)
        data.rename(columns={sample: 'DENOVO_' + sample}, inplace=True)
        column_names.append('^DENOVO_' + sample + '$')

    column_names += mq_column_names + COLUMNS_SET2
    data.rename(columns={'hits': 'Hits'}, inplace=True)

    if not all_data:  # get only experimental data
        data = data[data['Categories'].str.contains('Extra') & data['Extra'].notna()]

    # extract data from Leading razor protein (Extra column)
    data['HLA_affinity'] = data.apply(
        lambda x: 'SB' if x['HLA_rank'] < sb_threshold else ('WB' if x['HLA_rank'] < wb_threshold else None), axis=1)
    data['Gene_ID'] = data.apply(
        lambda x: set_from_extra(x['Sequence'], x['Extra'], x['Categories'], 'gene_id', unique=True),
        axis=1)
    data['Gene_Name'] = data.apply(
        lambda x: set_from_extra(x['Sequence'], x['Extra'], x['Categories'], 'gene_name', unique=True),
        axis=1)
    data['Pattern'] = data.apply(
        lambda x: set_from_extra(x['Sequence'], x['Extra'], x['Categories'], 'pattern', unique=True),
        axis=1)
    data['Out_Frame'] = data.apply(
        lambda x: set_from_extra(x['Sequence'], x['Extra'], x['Categories'], 'out_frame', unique=True),
        axis=1)
    data['DB'] = data.apply(
        lambda x: set_from_extra(x['Sequence'], x['Extra'], x['Categories'], 'db', unique=True),
        axis=1)
    data['Source_db'] = data.apply(
        lambda x: set_from_extra(x['Sequence'], x['Extra'], x['Categories'], 'source_db', unique=True),
        axis=1)
    data['Seq'] = data.apply(
        lambda x: multi_matching(x['Sequence'], x['Extra'], unique=True),
        axis=1)
    data['Chimera'] = data.apply(
        lambda x: 'y' if pd.notna(x['Seq']) and re.search(r'[A-Z][a-z]', x['Seq']) else 'n', axis=1)

    column_set = []
    for name in column_names:
        column_set += list(data.columns[data.columns.str.contains(name)])

    return data[column_set]

# end of  set_samples()


def main():
    parser = argparse.ArgumentParser(description="A script to filter experimental data")

    parser.add_argument('-i', required=True, help="input CSV file (the output file of integration_filter.py script)")
    parser.add_argument('-o', default='filtered_file4experimental.csv', required=False, help='output CSV file')
    parser.add_argument('-all', action='store_true', help='get all data, not only experimental (`Extra`)')
    parser.add_argument('-sb', default=0.5, required=False, help='strong binder (SB) threshold (default < 0.5)')
    parser.add_argument('-wb', default=2.0, required=False, help='weak binder (WB) threshold (default < 2.0)')

    args = parser.parse_args()

    input_file = args.i
    output_file = args.o
    all_data = args.all
    sb_threshold = args.sb
    wb_threshold = args.wb

    try:
        data = make_table(input_file, sb_threshold, wb_threshold, all_data)
        data.to_csv(output_file, sep='\t', index=False)
    except Exception as err:
        print("Something went wrong: {}".format(err))

    print('Extended integration filter: done')

# end of main()


if __name__ == '__main__':
    main()
