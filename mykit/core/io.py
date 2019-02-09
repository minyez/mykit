# -*- coding: utf-8 -*-
'''
defined classes for file manipulation
'''

import os
from shutil import copy2

# ==================== metaclasses ====================
class __file:
    '''Meta-metaclass for file reading and writing
    '''
    path = ''
    __path = ''
    __lines = []
    __fExist = False

    def __init(self, filePath):
        self.path = filePath
        self.__lines = []
        self.__fExist = False
    
    def __open(self):
        '''Open existing file
        '''
        if os.path.isfile(self.path):
            with open(self.__path, 'r') as hIn:
                self.__lines = hIn.readlines()
                self.__fExist = True
        else:
            self.__fExist = False


class file_in(__file):
    '''metaclass to create and manipulate input file of calculations

    Args:
        filePath (str): The name of the input file
    '''

    def __init__(self, filePath):
        self.__init(filePath)
        self.__open()
        self.iLines = self.__lines #: the lines of input file to write out

    def write(self, outPath=None, overWrite=False):
        '''Write to outPath

        Args:
            outPath (str) : the path of file to write the :data:`iLines`. Default use :data:`path`
                attribute
            overWrite (bool) : flag to overwrite if the file exists
        '''
        if outPath:
            __outpath = outPath
        else:
            __outpath = self.path

        if os.path.isfile(__outpath) and overWrite:
            pass
        else:
            copy2(__outpath, __outpath + '_bak')

        with open(__outpath, 'w') as h:
            for line in self.__lines:
                h.write(line.strip('\n') + '\n')


class file_out(__file):
    '''metaclass for output file of calculations

    Args:
        filePath (str): The name of the output file
    '''
    def __init__(self, filePath):
        self.__init(filePath)
        self.__open()
        self.oLines = self.__lines


# ========== Functions =============
def read_out(filePath):
    '''Read the output file 
    '''
    return file_out(filePath)

def set_in(filePath):
    '''Initialize the input file
    '''
    return file_in(filePath)
