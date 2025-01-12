"""CSVtools.py: The tools for working with CSV files"""

__author__ = "Dmitry Malko"


import re
import csv
import gzip


class CSV:
    @staticmethod
    def get_delimiter(file_path):
        sniffer = csv.Sniffer()
        if re.search(r'.gz$', file_path):
            with gzip.open(file_path, "rt") as f:
                header = f.readline()
        else:
            with open(file_path, "rt") as f:
                header = f.readline()

        delimiter = sniffer.sniff(header).delimiter
        if not re.search(r'[\t, ]', header):
            delimiter = '\t'  # fix the delimiter for one column file

        if delimiter not in ['\t', ',']:
            file_name = re.sub(r'.*/', '', file_path)
            raise ValueError('Can not determine CSV delimiter for {} file'.format(file_name))

        return delimiter

    # end of get_delimiter()

# end of class CSV
