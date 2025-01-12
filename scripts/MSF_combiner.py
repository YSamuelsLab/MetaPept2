#!/usr/bin/env python3

"""MSF_combiner.py: A script to organize FragPipe outputs"""

__author__ = "Dmitry Malko"


import re
import argparse
from glob import glob
import pandas as pd
import xml.etree.ElementTree as ET


class FragPipeCombiner:
    # INPUT
    _main_file = 'combined_peptide.tsv'
    _psm_file = 'psm.tsv'
    _xml_file = 'interact-*.pep.xml'  # for future projects
    _ion_file = 'ion.tsv'  # for future projects

    _column2rename_peptide = {
        'Peptide Sequence': 'Sequence',
        'Peptide Length': 'Length'
    }

    # DON'T PUT `Intensity` before `MaxLFQ Intensity` !!!
    _column2reorder_peptide = ['Spectral Count', 'MaxLFQ Intensity', 'Intensity', 'Match Type']

    _column2rename_msms = {
        'Peptide': 'Sequence',
        'Peptide Length': 'Length',
        'Calculated Peptide Mass': 'Mass',
        'Retention': 'Retention time',
        'PeptideProphet Probability': 'PeptideProphet'
    }

    # _dummy_column_peptide = ['Reverse', 'Leading razor protein', 'Potential contaminant']
    _dummy_column_peptide = []
    # _dummy_column_msms = ['Matches', 'Intensities']
    _dummy_column_msms = []

    # OUTPUT
    _output_peptide_file = 'peptides.txt'
    _output_msms_file = 'msms.txt'

    # XML (optional feature)
    _summary = 'interact_summary'
    _query = 'spectrum_query'
    _hit = 'search_hit'

    def __init__(self, input_dir):
        self._input = re.sub(r'/$', '', input_dir)
        self._pept_data = None
        self._msms_data = None

    def make_pep(self, input_file):
        main_file_name = '/'.join([self._input, input_file])
        self._pept_data = pd.read_csv(main_file_name, sep='\t', engine='python')

        if len(self._column2rename_peptide.keys()):
            self._pept_data = self._pept_data.rename(columns=self._column2rename_peptide)

        for col in self._column2reorder_peptide:  # columns reordering for MSFragger
            self._pept_data = self._pept_data.rename(columns=lambda x: re.sub(r'^(.*) ({})$'.format(col), r'\2 \1', x))

        if 'Intensity' not in self._pept_data.columns:  # SUM of Intensities for MSFragger
            columns = [x for x in self._pept_data.columns if re.match('Intensity', x)]
            self._pept_data['Intensity'] = self._pept_data[columns].sum(axis=1)

        for dummy_column in self._dummy_column_peptide:
            self._pept_data[dummy_column] = '.'

        return self._pept_data

    def make_msms(self):
        psm_file_name = '/'.join([self._input, '**', self._psm_file])

        msms_data = []
        for psm_file in glob(psm_file_name, recursive=True):
            data_item = pd.read_csv(psm_file, sep='\t', engine='python')

            # fit the data to MaxQuant output
            data_item['Spectrum File'] = data_item['Spectrum File'].apply(lambda x: re.sub(r'.*[/\\]', '', x))
            data_item['Raw file'] = data_item['Spectrum'].apply(lambda x: (re.split(r'\.', x))[0])
            data_item['Scan number'] = data_item['Spectrum'].apply(lambda x: (re.split(r'\.', x))[1])
            data_item['Delta score'] = data_item.apply(lambda x: x['Hyperscore'] - x['Nextscore'], axis=1)
            # data_item['PEP'] = data_item['PeptideProphet Probability'].apply(lambda x: 1 - x)

            msms_data.append(data_item)

        self._msms_data = pd.concat(msms_data)

        if len(self._column2rename_msms.keys()):
            self._msms_data = self._msms_data.rename(columns=self._column2rename_msms)

        for dummy_column in self._dummy_column_msms:
            self._msms_data[dummy_column] = '.'

        # let's make the data consistent
        if self._pept_data is None:
            raise ValueError('ERROR: there is no data to check the consistency, run make_pep() method before')
        n_before = len(self._msms_data['Sequence'].unique())
        self._msms_data = pd.merge(self._msms_data, self._pept_data['Sequence'], on=['Sequence'], how='inner')
        n_after = len(self._msms_data['Sequence'].unique())
        if n_before - n_after > 0:
            print('{} unique peptides were removed from psm.tsv data'.format(n_before - n_after))

        return self._msms_data

    def make_msms_xml(self):  # a method to parse XML output file
        psm_file_name = '/'.join([self._input, '**', self._psm_file])
        file_name = re.sub(r'.*/', '', psm_file_name)
        msms_data = []
        for psm_file in glob(psm_file_name, recursive=True):
            dir_path = re.sub(r'(.*)/.*', r'\1', psm_file)
            data_item = pd.read_csv(psm_file, sep='\t', engine='python')
            data_item['Spectrum File'] = data_item['Spectrum File'].apply(lambda x: re.sub(r'.*[/\\]', '', x))
            xml_files = list(data_item['Spectrum File'].unique())
            if len(xml_files) != 1:
                raise ValueError('ERROR: the file {} has no link to XML Spectrum File or multiple links'.format(psm_file))

            xml_file = xml_files[0]
            data_xml = ET.parse('/'.join([dir_path, xml_file]))
            xml_root = data_xml.getroot()
            # let's take into account XML NameSpace
            s = re.search(r'{([^{}]+)}', xml_root.tag)
            spectrum_query, xml_namespace = \
                ('.//ns:' + self._query, {'ns': s.group(1)}) if s else ('.//' + self._query, {})

            xml_data = []
            for query in xml_root.iterfind(spectrum_query, xml_namespace):
                spec = query.attrib['spectrumNativeID']
                s = re.search(r'\sscan=(\d+)', spec)
                scan = s.group(1) if s else 0
                pseq = ''
                search_hit = './/ns:' + self._hit if len(xml_namespace) else './/' + self._hit
                for hit in query.iterfind(search_hit, xml_namespace):
                    pseq = hit.attrib['peptide']
                    break

                xml_data.append({'Spectrum File': xml_file,'Peptide': pseq, 'Scan': scan})

            xml_data = pd.DataFrame(xml_data, index=None)
            data_item = pd.merge(data_item, xml_data, how='left', on=['Peptide'])
            msms_data.append(data_item)

        self._msms_data = pd.concat(msms_data)

        return self._msms_data

    def save(self):
        peptide_file_name = '/'.join([self._input, self._output_peptide_file])
        msms_file_name = '/'.join([self._input, self._output_msms_file])

        if len(self._pept_data.index):
            self._pept_data.to_csv(peptide_file_name, sep='\t', index=False)

        if len(self._msms_data.index):
            self._msms_data.to_csv(msms_file_name, sep='\t', index=False)

        return peptide_file_name, msms_file_name


def main():
    parser = argparse.ArgumentParser(description='A script to organize FragPipe outputs')
    parser.add_argument('-d', required=True, help='the input directory')
    parser.add_argument('-c', required=True, help='the name of combined_peptide.tsv file after contamination cleaning')

    args = parser.parse_args()
    input_dir = args.d
    input_file = args.c

    try:
        fp = FragPipeCombiner(input_dir)
        peptide_data = fp.make_pep(input_file)
        msms_data = fp.make_msms()
        peptide_file, msms_file = fp.save()
    except Exception as err:
        print('Something went wrong: {}'.format(err))
    else:
        print('{}, {} files are created'.format(peptide_file, msms_file))
        print('FragPipe file processing: done')

    return None

# end of main()


if __name__ == '__main__':
    main()
