# coding = utf-8
'''Class that manipulates WIEN2k main input struct file 
'''
from mykit.core.cell import Cell
from mykit.core.symmetry import get_spacegroup, get_sym_operations, standardize
from mykit.core.utils import get_cwd_name



class StructError(Exception):
    pass


class Struct(Cell):
    '''Class to manipulate WIEN2k master input case.struct file

    Args:
        latt, atoms, pos, kwargs: see docstring of Cell
        rmt (dict)
        r0 (dict)
    '''

    def __init__(self, latt, atoms, pos, rmt=None, r0=None, **kwargs):
        super(Struct, self).__init__(latt, atoms, pos, **kwargs)
        self._set_rmt_r0(rmt, r0)
    
    def _set_rmt_r0(self, rmt, r0):
        pass
    
    def __str__(self):
        return ''

    @classmethod
    def read_from_file(self, pathStruct=None):
        '''Read from an existing struct file

        Args:
            pathStruct (str): path to the file to read as WIEN2k struct
        '''
        p = pathStruct
        if p is None:
            p = get_cwd_name() + '.struct'
        
        with open(p, 'r') as h:
            slines = h.readlines()
        
        
        
        raise NotImplementedError
