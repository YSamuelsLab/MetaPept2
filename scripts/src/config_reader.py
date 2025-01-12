# -*- coding: utf-8 -*-
"""
Created on Thu June 15 10:27 2023

Usage: 	read config ('.ini') file and store arguments for the crIMP pipeline

Details: 

@author: Bracha Erlanger Avigdor

"""


import configparser
import sys


class ConfigReader:
    def __init__(self, file_path):
        self.file_path = file_path
        self.config = configparser.ConfigParser()
        self.config.read(self.file_path)

    def get_sections(self):
        return self.config.sections()

    def get_options(self, section):
        return self.config.options(section)

    def get_value(self, section, option):
        return self.config.get(section, option)

