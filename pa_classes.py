#!/usr/bin/env python
# coding=utf-8

# ====================================================
#
#     File Name : pa_classes.py
# Creation Date : 09-04-2018
#    Created By : Min-Ye Zhang
#       Contact : stevezhang@pku.edu.cn
#
# ====================================================

from __future__ import print_function, absolute_import
from pc_elements import common_read_chemical_formula
from pc_utils import common_run_calc_cmd
import sys, os
import subprocess as sp
from fnmatch import fnmatch

class abinit_input_files():

    '''
    Class for manipulating the abinit input file.
    '''


    def __init__(self, casename, formula, xc_type='pbe', pp_type='paw', **kwargs):
        '''
        Initialize the master input file i.e. casename.in file, along with 
        the input control casename.files
        '''

        self.pp_type_avail = {
                         'paw'  : 'xml'  ,
                         'ncpp' : 'psp8' ,
                         'usp'  : 'usp'  ,
                        }

        self.casename = casename
        self.abinit_cmd = None
        self.abinit_path = None

        self.__write_atom_info(formula)
        self.__write_pp_info(xc_type, pp_type)
        self.__set_files()



    def set_abinit_run(self, abinit_address_in=None, nproc=1, mpitype='mpiexec'):
        '''
        Initialize the running script of abinit calculation
        '''
        try:
            assert type(nproc)==int
            assert nproc > 0
        except AssertionError:
            raise TypeError("nproc should be a positive integer.")

        if abinit_address_in is None or not os.path.exists(abinit_address_in):
            if not os.path.exists(abinit_address_in):
                print("WARNING: Specified executive not found: %s\nSearching from PATH..." % abinit_address_in)
            try:
                abinit_address = sp.check_output(['which','abinit']).strip()
            except sp.CalledProcessError:
                raise ValueError("ABINIT executive not found. Please figure it")
            else:
                print("Found ABINIT executive at %s" % abinit_address)
                self.abinit_path = abinit_address
        elif os.path.exists(abinit_address_in):
            self.abinit_path = os.path.abspath(abinit_address_in)

        self.nproc = nproc

        if nproc==1:
            self.abinit_cmd = self.abinit_path + ' < ' + self.controlfiles
        else:
            self.abinit_cmd = mpitype + ' -np ' + str(nproc) + ' ' + self.abinit_path + \
                              ' < ' + self.controlfiles


    def run_abinit(self, fout=None, ferr=None):
        '''
        Perform calculation after setting up abinit_cmd
        '''
        try:
            assert self.abinit_cmd is not None
        except AssertionError:
            raise ValueError('The abinit_cmd attribute has not been initialized yet.')
        else:
            common_run_calc_cmd(self.abinit_cmd, fout, ferr)


    def __write_atom_info(self, formula):
        '''
        Parameters:
            formula: (str)
                the non-reduced chemical formula for the system to calculate. 
        '''
        self.atom_type, self.natom = common_read_chemical_formula(formula)
        self.natoms = sum(self.natom)


    def __write_pp_info(self, xc_type, pp_type):
        '''
        TODO: support more functionals.
        '''

        if xc_type.lower() == 'pbe':
            self.xc_type = 'pbe'
        else:
            # Functionals other than PBE is not supported currently
            raise ValueError("Functional other than PBE not supported yet.")

        # get the environment variable for directores which stores the psp8 NCPP
        # and the PAW dataset
        # PAW should be a version of JTH PAW dataset, while NCPP is the nc-sr_pbe_standard_psp8
        try:
            self.abinit_pp_path_list = os.environ['ABINIT_PP_PATH'].split(':')
        except KeyError:
            raise KeyError('ABINIT_PP_PATH should be set in the environment variable.')
        else:
            if len(self.abinit_pp_path_list) == 1 and len(self.abinit_pp_path_list[0]) == 0:
                raise ValueError('ABINIT_PP_PATH is probably not correctly defined. Please check.')

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
            raise ValueError("Corresponding directory of PP type is not found in environment")

        # find the PP for the atoms. For NCPP and JTH-PAW, the naming convention are both "Si.xxxxx"
        self.atom_pp_list = []
        for at in self.atom_type:
            for pp in os.listdir(current_pp_path):
                if fnmatch(pp, at+'.*'):
                    self.atom_pp_list.append(current_pp_path+'/'+pp)
                    break


    def __set_files(self):
        '''
        Set input file attributes and Generate input control file by the casename i.e. casename.files
        '''
        casename = self.casename
        self.infile = casename+'.in'
        self.outfile = casename+'.out'
        self.logfile = casename+'.log'
        self.controlfiles = casename + '.files'
        with open(self.controlfiles,'w') as h_files:
            h_files.write(self.infile+'\n')
            h_files.write(self.outfile+'\n')
            h_files.write(casename+'i\n')
            h_files.write(casename+'o\n')
            h_files.write(casename+'\n')
            for pp in self.atom_pp_list:
                h_files.write(pp+'\n')



