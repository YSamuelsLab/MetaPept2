#!/usr/bin/env python3

import argparse
import re
import pandas as pd
from CSVtools import CSV

# DEFAULT FILTERING VALUES

# combined hits
# ALC_combined_sufficient = 80  # no other filters for this threshold
# ALC_combined = 70  # the additional filters for this threshold are below
# Coverage = 80
# Delta = 10

# unique hits
# ALC_denovo = 80
# Q_denovo = 0.1
# HLA_rank = 2.0


def get_best_coverage(row):
    best = row.filter(regex=r'^coverage').max()
    return best

# end of get_best_coverage()


def get_best_delta(row):
    best = row.filter(regex=r'^Delta_score').max()
    return best

# end of get_best_delta()


def combined_filter(input_file, ALC_combined_sufficient, ALC_combined, Coverage, Delta, Hyperscore, Deltascore):
    data = pd.read_csv(input_file, sep=CSV.get_delimiter(input_file), engine='python')
    data['Best_coverage'] = data.apply(get_best_coverage, axis=1)
    data['Best_delta_score'] = data.apply(get_best_delta, axis=1)

    if True in data.columns.str.contains('Hyperscore'):
        # MSFragger data
        data = data[(data['Best_ALC'] >= ALC_combined_sufficient) | (
                (data['Best_ALC'] < ALC_combined_sufficient) &
                (data['Best_ALC'] >= ALC_combined) &
                (data['BestHit_Hyperscore'] >= Hyperscore) &
                (data['BestHit_Deltascore'] >= Deltascore))
        ]
    else:
        # MaxQuant data
        data = data[(data['Best_ALC'] >= ALC_combined_sufficient) | (
            (data['Best_ALC'] < ALC_combined_sufficient) &
            (data['Best_ALC'] >= ALC_combined) &
            (data['Best_coverage'] >= Coverage) &
            (data['Best_delta_score'] >= Delta))
        ]

    data = data.rename(columns={'HLA_Allele': 'HLA_allele'})
    data = data.rename(columns=lambda x: re.sub(r'(HLA_.*)_EL_Rank', r'\1', x))  # for human
    data = data.rename(columns=lambda x: re.sub(r'(H_2_.*)_EL_Rank', r'\1', x))  # for mouse
    data.drop(columns=['Best_PRISM_Replica', 'Filtered_HLA_allele'], inplace=True)
    data['HLA_allele'] = data['HLA_allele'].apply(lambda x: re.sub(r':', '', re.sub(r'-', '_', x)))
    data = data.rename(columns=lambda x: re.sub(r'(.*)_HLA_.*', r'\1', x))  # for human
    data = data.rename(columns=lambda x: re.sub(r'(.*)_MHC_.*', r'\1', x))  # for mouse

    data['Integration'] = 'Combined'

    return data

# end of combined_filter()


def denovo_filter(input_file, ALC_denovo, Q_denovo, HLA_rank):
    data = pd.read_csv(input_file, sep=CSV.get_delimiter(input_file), engine='python')
    data = data.loc[(data['Best_ALC'] >= ALC_denovo) & (data['Best_Q'] <= Q_denovo)]
    data = data.loc[data['netMHC_rank'] < HLA_rank]

    data = data.rename(columns={'netMHC_rank': 'HLA_rank'})
    data = data.rename(columns={'Intensity_Sum': 'Intensity'})
    data = data.rename(columns=lambda x: re.sub(r'(Intensity_[^_]+_)[^0-9]+(\d+)', r'\1\2', x))
    data.drop(columns=['Location_count', 'Genome', 'Top_location_count_no_decoy', 'Filtered_HLA_allele'], inplace=True)
    data['Integration'] = 'Denovo'

    return data

# end of denovo_filter()


def msf_filter(input_file, hyperscore, deltascore):
    if input_file is None:
        return pd.DataFrame()

    data = pd.read_csv(input_file, sep=CSV.get_delimiter(input_file), engine='python')
    data['CDS'] = data['CDS'].astype(str)
    data['nuORFs'] = data['nuORFs'].astype(str)

    data = data.loc[(data['BestHit_Hyperscore'] >= hyperscore) & (data['BestHit_Deltascore'] >= deltascore)]

    data = data.rename(columns={'HLA_Allele': 'HLA_allele'})
    data = data.rename(columns=lambda x: re.sub(r'(HLA_.*)_EL_Rank', r'\1', x))  # for human
    data = data.rename(columns=lambda x: re.sub(r'(H_2_.*)_EL_Rank', r'\1', x))  # for mouse
    data['HLA_allele'] = data['HLA_allele'].apply(lambda x: re.sub(r':', '', re.sub(r'-', '_', x)))
    data = data.rename(columns=lambda x: re.sub(r'(.*)_HLA_.*', r'\1', x))  # for human
    data = data.rename(columns=lambda x: re.sub(r'(.*)_MHC_.*', r'\1', x))  # for mouse

    data['Categories'] = data.apply(lambda x: 'CDS' if x['CDS'] != 'nan' else 'nuORFs' if x['nuORFs'] != 'nan' else 'Extra', axis=1)
    data['Integration'] = 'MSFragger'

    return data

# end of denovo_filter()


def make_output(output_file, com_data, denovo_data, imp_data):
    data = pd.concat([com_data, denovo_data, imp_data], axis=0, ignore_index=True)
    data = data.drop(columns=['Source_File'])

    # set `Sequence` as the first column
    col_seq = data.pop('Sequence')

    data.insert(0, col_seq.name, col_seq)
    data.to_csv(output_file, sep='\t', index=False)
    return data

# end of make_output()


def main():
    parser = argparse.ArgumentParser(description="A script to filter integrated data based on Q, ALC, Coverage and Delta")

    parser.add_argument('-com', required=True, help="combined data file")
    parser.add_argument('-denovo', required=True, help="denovo unique data file")
    parser.add_argument('-imp', default=None, help="IMP unique data file")
    parser.add_argument('-o', default='filtered_file.csv', required=False, help='output file with filtered data')
    parser.add_argument('-alc_suff', default=80, required=False, help='sufficient ALC threshold for combined hits')
    parser.add_argument('-alc_comb', default=70, required=False, help='ALC threshold for combined hits')
    parser.add_argument('-cov_comb', default=80, required=False, help='Coverage threshold for combined hits')
    parser.add_argument('-delta_comb', default=10, required=False, help='Delta score threshold for combined hits')
    parser.add_argument('-alc_denovo', default=80, required=False, help='ALC threshold for denovo unique hits')
    parser.add_argument('-q_denovo', default=0.1, required=False, help='Q (FDR) threshold for denovo unique hits')
    parser.add_argument('-rank_denovo', default=2.0, required=False, help='Rank threshold for denovo unique hits')
    parser.add_argument('-hyper_msf', default=20, required=False, help='Hyperscore threshold for MSFragger unique hits')
    parser.add_argument('-delta_msf', default=4, required=False, help='Deltascore threshold for MSFragger unique hits')

    args = parser.parse_args()

    combined_file = args.com
    denovo_file = args.denovo
    imp_file = args.imp
    output_file = args.o
    alc_suff = args.alc_suff
    alc_comb = args.alc_comb
    cov_comb = args.cov_comb
    delta_comb = args.delta_comb
    alc_denovo = args.alc_denovo
    q_denovo = args.q_denovo
    rank_denovo = args.rank_denovo

    hyper_msf = args.hyper_msf
    delta_msf = args.delta_msf

    try:
        denovo_data = denovo_filter(denovo_file, alc_denovo, q_denovo, rank_denovo)
        com_data = combined_filter(combined_file, alc_suff, alc_comb, cov_comb, delta_comb, hyper_msf, delta_msf)
        imp_data = msf_filter(imp_file, hyper_msf, delta_msf)
        make_output(output_file, com_data, denovo_data, imp_data)
    except Exception as err:
        print("Something went wrong: {}".format(err))

    print('Integration filter: done')

# end of main()


if __name__ == '__main__':
    main()
