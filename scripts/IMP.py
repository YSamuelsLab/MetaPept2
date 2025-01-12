#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 11:30:27 2023


Master script running IMP program.

Usage: 
    

Details: 
    
    
@author: Bracha Erlanger Avigdor

"""

__version__ = '1.0.1'

import configparser
import pathlib
import argparse
import sys
import os
import numpy as np
import pandas as pd

sys.path.insert(0, './src')

from src.peptides import peptides
from src.NetMHCpan import netMHCpan
from src.msms import msms
from src.merge_tables import mergeTables
from src.filter_tables import filterTables


def is_valid_path(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return arg  # return arg as is


def main(args):
    """
    Main program to call the IMP pipeline
    """

    print('calling peptide')
    # read peptide data and add I to L data and perform database search
    pp = peptides(args)

    print('calling netMHCpan')
    # call NetMHCpan
    nm = netMHCpan(args)
    nm.parse_netMHCpan()

    print('calling msms')
    # read MS-MS data
    ms = msms(args)

    print('merging tables')
    mt = mergeTables(args)
    mt.merge_MQ_netMHCpan_tables()

    print('filtering')
    ft = filterTables(args)
    ft.filter_MQ_netMHCpan_peptides()
    exit()


def make_parser():
    # define lambda function to validate file/path exists to use in add_argument
    is_valid = type = lambda x: is_valid_path(parser, x)

    parser = argparse.ArgumentParser(description="Immunopeptidomics pipeline",
                                     epilog="Samuels Lab 2023", prog='IMPly')

    _input = parser.add_argument_group('input options')

    _input.add_argument('-i', '--input_folder', metavar='', type=is_valid, default="./input",
                        help="folder location of input file from MaxQuant (default: %(default)s)")
    _input.add_argument('-a', '--alleles', metavar='', type=str, default="HLA.txt",
                        help="HLA alleles used in the experiment (default: %(default)s)")
    _input.add_argument('-p', '--peptides_file', metavar='', type=str, default="peptides.txt",
                        help="peptides file from MaxQuant (default: %(default)s)")
    _input.add_argument('-m', '--msms_file', metavar=',', type=str, default="msms.txt",
                        help="MS/MSfile from MaxQuant (default: %(default)s)")
    _input.add_argument('-s', '--sample_desc_file', metavar='Sample description file', type=str,
                        default="sample_description.csv",
                        help="sample description file (default: %(default)s) ")

    _input.add_argument('--dummy', action='store_true',
                        help="DO NOT predict binding")
    _input.add_argument('--max_len', metavar='maximum peptide length', type=int, default=14,
                        help="maximum peptide length")

    _input.add_argument('-d', '--database_folder', metavar='', type=is_valid, default="./db",
                        help="folder location of database files to search (default: %(default)s) ")
    _input.add_argument('-c', '--CDS_fasta_file', metavar='', type=str, default="CDS.fasta",
                        help="CDS database file in fasta format (default: %(default)s) ")
    _input.add_argument('-n', '--nuORFdb_fasta_file', metavar='', type=str, default="nuORFdb.fasta",
                        help="nuORF database file in fasta format (default: %(default)s) ")

    _output = parser.add_argument_group('output options')

    _output.add_argument('-o', '--output_folder', metavar='', type=is_valid, default="./output",
                         help="folder location for output files (default: %(default)s)")

    _netMHCpan = parser.add_argument_group('netMHCpan options')

    _netMHCpan.add_argument('-b', '--netMHCpan_path', metavar='', type=is_valid,
                            default="./tools/netMHCpan-4.1/Linux_x86_64/",
                            help="netMHCpan bin, path (default: %(default)s)")
    _netMHCpan.add_argument('-t', '--tmpdir', metavar='', type=is_valid, default="/tmp",
                            help="writable temp directory path (default: %(default)s)")

    parser.add_argument('-v', '--version', action='version', version="v%s" % (__version__))

    return parser


if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args()

    # call main
    main(args)
    print('IMP: done')
