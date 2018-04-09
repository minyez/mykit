#!/usr/bin/env python
# coding=utf-8

# ====================================================
#
#     File Name : pa_input.py
# Creation Date : 09-04-2018
#    Created By : Min-Ye Zhang
#       Contact : stevezhang@pku.edu.cn
#
# ====================================================

from __future__ import print_function, absolute_import
from pc_elements import common_read_chemical_formula
from pc_elements import Periodic_Table as PT
import sys, os
from fnmatch import fnmatch

class abinit_input_files():

    '''
    Class for manipulating the abinit input file.
    '''

    pp_type_avail = {
                     'paw' : 'xml',
                     'ncpp': 'psp8'
                    }

    # get the environment variable for directores which stores the psp8 NCPP
    # and the PAW dataset
    # PAW should be a version of JTH PAW dataset, while NCPP is the nc-sr_pbe_standard_psp8
    abinit_pp_path_list = os.environ['ABINIT_PP_PATH'].split(':')


    def __init__(self, casename, formula, xc_type='pbe', pp_type='paw'):
        '''
        initialize the master input file i.e. casename.in file, along with 
        the input files control casename.files
        '''
        self.casename = casename
        self.__write_atom_info(formula)
        self.__write_pp_info(xc_type, pp_type)
        self.__set_files()


    def __write_atom_info(self, formula):
        '''
        Parameters:
            formula: string
                the non-reduced chemical formula for the system to calculate. 
        '''
        self.atom_type, self.natom = common_read_chemical_formula(formula)
        self.natoms = sum(self.natom)

 


    def __write_pp_info(self, xc_type, pp_type):
        '''
        TODO: only support PBE.
        '''

        if xc_type.lower() == 'pbe':
            self.xc_type = 'pbe'
        else:
            # Functionals other than PBE is not supported currently
            raise ValueError("Functional other than PBE not supported yet.")

        # set pseudopotential path, both ncpp and paw

        for pp_path in self.abinit_pp_path_list[::-1]:
            if fnmatch(pp_path, '*JTH*'):
                self.abinit_paw_path = pp_path
            if fnmatch(pp_path, '*pbe_sr_s*'):
                self.abinit_ncpp_path = pp_path
            if fnmatch(pp_path, '*nc-sr_pbe_standard_psp8*'):
                self.abinit_ncpp_path = pp_path

        if pp_type.lower() not in self.pp_type_avail:
            raise ValueError("Unavailable pseudopotential type.")
        elif pp_type.lower() == 'paw' and self.abinit_paw_path is not '':
            current_pp_path = self.abinit_paw_path
        elif pp_type.lower() == 'ncpp' and self.abinit_ncpp_path is not '':
            current_pp_path = self.abinit_ncpp_path
        else:
            raise TypeError("Corresponding directory of PP type is not found in environment")

        # find the PP for the atoms. For NCPP and JTH-PAW, the naming convention are both "Si.xxxxx"
        self.atom_pp_list = []
        for at in self.atom_type:
            for pp in os.listdir(current_pp_path):
                if fnmatch(pp,at+'.*'):
                    self.atom_pp_list.append(current_pp_path+'/'+pp)
                    break


    def __set_files(self):
        casename = self.casename
        with open(casename+'.files','w') as h_files:
            h_files.write(casename+'.in\n')
            h_files.write(casename+'.out\n')
            h_files.write(casename+'i\n')
            h_files.write(casename+'o\n')
            h_files.write(casename+'\n')
            for pp in self.atom_pp_list:
                h_files.write(pp+'\n')
