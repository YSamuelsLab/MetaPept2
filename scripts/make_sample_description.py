#!/usr/bin/env python3

import re
import argparse
import glob
import pandas as pd


def description(input_dir, output_file):
    source_files = set()
    file_names = glob.glob(re.sub(r'/$', '', input_dir) + '/**/*.csv.gz', recursive=True)
    if len(file_names):
        for full_name in file_names:
            name = re.sub(r'.*/', '', full_name)
            if re.search(r'\.csv\.gz.*\.csv\.gz$', name) or re.search(r'^i', name):
                # it is not a PRISM input file
                continue

            print('file processing: {} '.format(name), end='')
            data = pd.read_csv(full_name, compression='gzip', sep=',', quotechar='"', low_memory=False)
            source_files.update(data['Source File'])
            print('...OK')
    else:
        file_names = glob.glob(re.sub(r'/$', '', input_dir) + '/**/*de novo*.csv', recursive=True)
        for full_name in file_names:
            name = re.sub(r'.*/(.*/)', r'\1', full_name)
            print('file processing: {} '.format(name), end='')
            data = pd.read_csv(full_name, sep=',', quotechar='"', low_memory=False)
            source_files.update(data['Source File'])
            print('...OK')

    desc = pd.DataFrame(sorted(source_files), columns=['Source_File'])
    desc['Sample_Name'] = ''
    desc['Sample_Replica'] = ''
    desc['Sample_Type'] = ''
    desc['Group'] = ''
    desc['Experiment'] = ''

    desc.to_csv(output_file, sep='\t', index=False)

    return desc['Source_File'].to_list


def main():
    parser = argparse.ArgumentParser(description="A script for building a sample description file "
                                                 "for the PRIMS pipeline")

    parser.add_argument('-i', required=True, help="input directory with `*.csv.gz` or `*de novo*.csv` files "
                                                  "(they will be found recursively)")
    parser.add_argument('-o', default='sample_description.csv', required=False, help='the sample description file')

    args = parser.parse_args()

    in_dir = args.i
    out_file = args.o

    description(in_dir, out_file)

    print('...done')


if __name__ == '__main__':
    main()
