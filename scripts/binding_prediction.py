#!/usr/bin/env python3

import os
import re
import argparse
import pandas as pd
from glob import glob
import shutil
from typing import Callable


list_of_tools = ['Dummy', 'netMHCpan', 'netMHCIIpan']


class ToolNotFoundError(ValueError):
    def __init__(self):
        print('ERROR: unrecognizable tool!')


class FactorySubject(object):  # Factory pattern to manage classes
    columns = ['Sequence', 'HLA allele', 'netMHC % rank']

    def __init__(self):
        self.alleles = []

    def set(self, alleles=None):
        self.alleles = alleles if isinstance(alleles, list) else self.alleles


class Dummy(FactorySubject):
    # Dummy predictor to make PRISM output compatible with MetaPept pipeline in the 'no binding prediction' mode
    allele = 'HLA-A*00:00'

    def __init__(self):
        super().__init__()

    def predict(self, pep_list):
        my_data = pd.DataFrame(pep_list, columns=['Sequence'])
        my_data['HLA allele'] = self.allele
        my_data['netMHC % rank'] = 0
        my_data[self.allele] = 0

        return my_data


class netMHCpan(FactorySubject):
    path2tool = '../tools/netMHCpan/netMHCpan'

    def __init__(self):
        super().__init__()

    def predict(self, pep_list):
        # TODO: write code for running netMHCpan
        ...


class netMHCIIpan(FactorySubject):
    path2tool = '../tools/netMHCIIpan/netMHCIIpan'

    def __init__(self):
        super().__init__()

    def predict(self, pep_list):
        # TODO: write code for running netMHCIIpan
        ...


class ToolFactory(object):
    @staticmethod
    def get(class_name: str) -> object:
        if type(class_name) is not str:
            raise ValueError("class_name must be a string!")

        raw_subclasses_ = FactorySubject.__subclasses__()
        classes: dict[str, Callable[..., object]] = {c.__name__:c for c in raw_subclasses_}
        my_class = classes.get(class_name, None)
        if my_class is not None:
            return my_class

        raise ToolNotFoundError()


class BindingPredictor:
    def __init__(self, dir_path):
        self._files = []
        for file in glob(re.sub(r'/$', '', dir_path) + '/**/*.pep.annotated.csv.gz', recursive=True):
            if os.path.isfile(file):
                shutil.copy2(file, file + '.bak')
                self._files.append(file)

    #@staticmethod
    def run(self, tool_name, allele_file):
        tool_class = ToolFactory.get(tool_name)
        tool = tool_class()

        if allele_file:
            alleles = pd.read_csv(allele_file, header=None)
            alleles = alleles.iloc[:, 0].values
            tool.set(alleles)

        peptides = set()  # the full set of unique peptides
        for file in self._files:
            data = pd.read_csv(file, compression='gzip', sep=',', header=0)
            pept = list(data['Sequence'])
            peptides.update(pept)

        prediction = tool.predict(peptides)

        upd_files = []
        for file in self._files:
            data = pd.read_csv(file, compression='gzip', sep=',', header=0)
            data = pd.merge(data, prediction, on=['Sequence'], how='left')
            data.to_csv(file, compression='gzip', sep=',')

            print('{} is updated'.format(re.sub(r'.*/', '', file)))
            upd_files.append(file)

        return upd_files


def main():
    parser = argparse.ArgumentParser(description='A script for running MHC/peptide binding prediction tools '
                                     'on PRISM output files `*.pep.annotated.csv.gz`')
    parser.add_argument('-i', required=True,
                        help='the directory with PRISM output files `*.pep.annotated.csv.gz`')
    parser.add_argument('-t', choices=list_of_tools, default='Dummy',
                        help='a binding prediction tool; default: no binding prediction - the script will add pseudo '
                             'allele HLA-A*00:00 with prediction rank = 0 to make PRISM output compatible with '
                             'MetaPept pipeline')
    parser.add_argument('-a', default='', help='the file with MHC alleles')

    args = parser.parse_args()
    input_dir = args.i
    tool = args.t
    allele_file = args.a

    try:
        bp = BindingPredictor(input_dir)
        bp.run(tool, allele_file)
    except Exception as err:
        print("Something went wrong:", err)

    print('...done')

# end of main()


if __name__ == '__main__':
    main()
