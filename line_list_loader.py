from astropy.io import fits
import numpy as np
import pandas as pd
import os
import json


class LineListLoader:
    # Base class for specific line list loaders
    def __init__(self, filename):
        self.line_list_file = filename
        self.columns = pd.Index()

    def search_column_name(self, name):
        column_name_list = self.columns.array[self.columns.str.contains(name, case=False)]
        return column_name_list[0] if len(column_name_list) == 1 else None


class LineListLoader_Wise(LineListLoader):

    def __init__(self, filename='line_lists/wise/line_list_wise.txt'):
        self.line_list_file = filename
        self.dataframe = pd.read_csv(self.line_list_file, sep='\s\s+', engine='python')
        self.columns = self.dataframe.columns
        self.single_species = self.list_single_species()

    def list_single_species(self):
        '''Returns a DataFrame filtered by single species'''
        species_column_name = self.search_column_name('species')
        assert species_column_name is not None
        single_species_filter = self.dataframe[species_column_name].str.split(',').apply(len) == 1
        single_species = self.dataframe[single_species_filter]
        return single_species

    def list_VALD_depths(self, min_depth, max_depth):
        '''Returns a single_species DataFrame filtered by min and max VALD_Depths'''
        single_species = self.list_single_species()
        species_column_name = self.search_column_name('species')
        VALD_Depth_column = single_species[species_column_name].str.extract('\((0.\d+)\)', expand=False).astype(float)
        VALD_Depth_filter = (VALD_Depth_column >= min_depth) & (VALD_Depth_column <= max_depth)
        return single_species[VALD_Depth_filter]


class LineListLoader_VALD(LineListLoader):

    def __init__(self, filename='line_lists/vald/4300_6650/4300_6650.lin'):
        self.line_list_file = filename
        columns = ['Spec Ion', 'WL_vac(A)', 'Excit(eV)', 'Vmic', 'log gf*', 'Rad. damping', 'Stark damping', 'Waals damping', 'Lande factor', 'Central depth', 'References']
        self.dataframe = pd.read_csv(self.line_list_file, sep=',(?= |-|\d)', names=columns, usecols=range(11), skiprows=3, skipfooter=106, quoting="'", engine='python')
        self.dataframe['Spec Ion'] = self.dataframe['Spec Ion'].str.replace("'", '')
        self.dataframe['References'] = self.dataframe['References'].str.replace("'", '')
        self.columns = self.dataframe.columns

