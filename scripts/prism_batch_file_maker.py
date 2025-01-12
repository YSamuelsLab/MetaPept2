#!/usr/bin/env python3

import re
import os
import argparse
import gzip
import glob
import shutil
import pathlib

LSF_MODULES = ['R/4.1.2.rstudio-foss-2021b', 'jre/8.121']  # TODO: add a command line option to set arbitrary modules
PRISM_VERSION = 'tools/Peptide-PRISM/Prism_2023-01-16'
DE_NOVO_PEPTIDES = 'de novo peptides.csv'
ALL_CANDIDATES = 'all de novo candidates.csv'

# The hash below describes the default PRISM categories and aliases for them:
PRISM_CATS = {
    'frameshift': 'CDS,UTR5,OffFrame,UTR3,ncRNA,Frameshift,Intronic,Intergenic',
    'prio1': 'CDS,UTR5,OffFrame,UTR3,ncRNA,Frameshift,Intronic,Intergenic',  # the same as 'frameshift'
    'prio2': 'CDS,Extra,UTR5,OffFrame,UTR3,ncRNA,Intronic,Intergenic',
    'prio3': 'CDS,UTR5,OffFrame,UTR3,ncRNA,Extra,Intronic,Intergenic',
}


def set_categories(cats):
    aliases = set()
    for cat in cats:
        alias = re.sub(r'[^A-Za-z0-9_,=]', '_', cat)
        if re.search('=', alias):
            alias, category = re.split('=', alias)
            if alias in PRISM_CATS:
                print('WARNING: "{}" alias for the default category was rewritten!'.format(alias))
            PRISM_CATS[alias] = category
            aliases.add(alias)
        else:
            if alias in PRISM_CATS:
                aliases.add(alias)
            else:
                raise ValueError('"{}" alias has bad format'.format(cat))

    return sorted(aliases)


def set_netMHCpan(directory='./'):
    for filepath in pathlib.Path(directory).glob('*/**/netMHCpan'):
        if filepath.is_dir():
            continue
        
        abs_path = filepath.absolute()
        data = ''  # just to be sure that netMHCpan correctly installed
        with open(filepath, 'rt') as f:
            data = f.read()

        if not re.search(r'setenv\s+NMHOME\s+{}'.format(os.path.dirname(abs_path)), data):
            data = re.sub(r'setenv\s+NMHOME.*', 'setenv NMHOME {}\n'.format(os.path.dirname(abs_path)), data)
            with open(filepath, 'wt') as f:
                f.write(data)
            print('Warning: netMHCpan file was modified to set a correct pathway.')

        return abs_path

    return


def main():
    parser = argparse.ArgumentParser(description="PRISM batch runner")

    parser.add_argument('-i', metavar='dirname', required=True, help='input directory')
    parser.add_argument('-o', metavar='dirname', default='output', required=False, help='output directory')
    parser.add_argument('-r', metavar='filename', default='prism.run_batch.sh', required=False, help='batch file to run PRISM')
    parser.add_argument('-j', metavar='args', default='-Xmx60g -Xms6g', required=False,
                        help='java arguments (default "-Xmx60g -Xms6g")')
    parser.add_argument('-threads', metavar='n', required=False, help='The number of threads to use for computations')
    parser.add_argument('-g',  metavar='name', default='h.ens90', required=False, help='Genomic name')
    parser.add_argument('-extra', metavar='filename', required=False, help='Fasta file containing extra sequences')
    parser.add_argument('-hla', metavar='filename', required=False, help='File containing HLA alleles')
    parser.add_argument('-netmhc', metavar='path', required=False, help='Command to call netMHCpan')
    parser.add_argument('-cat', nargs='+', required=True,
                        help='list of aliases for PRISM categories (default aliases: '
                             'frameshift=CDS,UTR5,OffFrame,UTR3,ncRNA,Frameshift,Intronic,Intergenic, '
                             'prio1=CDS,UTR5,OffFrame,UTR3,ncRNA,Frameshift,Intronic,Intergenic, '
                             'prio2=CDS,Extra,UTR5,OffFrame,UTR3,ncRNA,Intronic,Intergenic, '
                             'prio3=CDS,UTR5,OffFrame,UTR3,ncRNA,Extra,Intronic,Intergenic)')
    parser.add_argument('-lsf', action='store_true', help='add modules to run on the LSF cluster')

    args = parser.parse_args()
    input_dir = args.i
    output_dir = args.o
    batch_file = args.r
    java_args = args.j
    lsf = args.lsf

    genome = args.g
    threads = args.threads  #
    extra_file = args.extra
    hla_file = args.hla
    netmhc = args.netmhc
    cats = args.cat

    cats = set_categories(cats)
    netMHCpan_path = set_netMHCpan()
    export_path = 'export PATH=$PATH:{}'.format(os.path.dirname(netMHCpan_path))
    if lsf:  # add additional module to run on the LSF cluster
        for lsf_module in LSF_MODULES:
            export_path += '\nmodule load ' + lsf_module

    prism = 'java {} -jar {}.jar'.format(java_args, PRISM_VERSION)
    prism += ' -nthreads {}'.format(threads) if threads else ''
    prism += ' -netmhc {}'.format(netmhc) if netmhc else ''
    prism += ' -g {}'.format(genome) if genome else ''
    # prism += ' -hla {}'.format(hla_file) if hla_file else ''  # HLA file path can be set directly in the batch file
    prism += ' -unidentified -deltaNextStats -D'

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    dir_list = os.listdir(input_dir)

    with open(batch_file, 'w') as run:
        print('Creating batch file...')
        print(export_path, file=run)

        for cat in cats:
            for name in dir_list:
                if os.path.isdir('/'.join([input_dir, name])):
                    new_hla_file = '/'.join([output_dir, cat + '.' + name]) + '.hla'
                    if hla_file:
                        shutil.copy(hla_file, new_hla_file)
                    else:
                        hla_files = glob.glob('/'.join([input_dir, name, '*.hla']))
                        if len(hla_files):
                            print('Found HLA file in {}'.format('/'.join([input_dir, name])))
                            shutil.copy(hla_files[0], new_hla_file)
                        else:
                            print("No HLA file!")
                            exit(1)

                    all_candidates_file = '/'.join([input_dir, name, ALL_CANDIDATES])
                    new_all_candidates_file = output_dir + '/' + cat + '.' + name + '.csv.gz'
                    with gzip.open(new_all_candidates_file, 'wb') as fo:
                        with open(all_candidates_file, 'rb') as fi:
                            fo.write(fi.read())

                    de_novo_pept_file = '/'.join([input_dir, name, DE_NOVO_PEPTIDES])
                    new_de_novo_pept_file = output_dir + '/i' + cat + '.' + name + '.csv.gz'
                    with gzip.open(new_de_novo_pept_file, 'wb') as fo:
                        with open(de_novo_pept_file, 'rb') as fi:
                            fo.write(fi.read())

                    prism_output = output_dir + '/' + cat + '.' + name + '.out'
                    prism_error = output_dir + '/' + cat + '.' + name + '.err'

                    category = ' -cat '
                    if cat in PRISM_CATS:
                        category += PRISM_CATS[cat]
                    else:
                        raise ValueError('"{}" alias has bad format'.format(cat))

                    extra = ''
                    if re.search(r'Extra', category):
                        if extra_file:
                            extra = ' -extra {}'.format(extra_file)
                        else:
                            print('Error: you use Extra category without extra file!')
                            exit(1)
                    command = prism + extra + category + ' -in ' + new_all_candidates_file + ' > ' + prism_output + ' 2> ' + prism_error
                    print(command, file=run)
                    print('{} -> {}'.format(name, cat))

    print('PRISM batch file maker: done')

# end of main()


if __name__ == '__main__':
    main()
