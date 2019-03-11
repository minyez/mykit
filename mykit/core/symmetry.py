# coding = utf-8
'''class and functions related to crystal symmetry
'''

import numpy as np
import spglib as spg

from mykit.core.cell import Cell
from mykit.core.numeric import prec


class SymmetryError(Exception):
    pass


class symmetry(prec):
    '''the class for symmetry information of crystal, powered by spglib

    Args:
        cell (Cell or its subclasses)
    '''
    
    def __init__(self, cell):

        try:
            assert isinstance(cell, Cell)
        except AssertionError:
            raise SymmetryError("The input should be instance of Cell or its subclass.")
        try:
            assert cell.unit == "D"
        except AssertionError:
            raise SymmetryError("The coordiante system should be direct. Cartisian found.")
            
        # convert to direct coordinate system
        self._cell = cell
        self.spaceGroup = spg.get_spacegroup(self._cell.get_spglib_input(), \
            symprec=self._symprec)
    
    def ibzkpt(self, kgrid, shift=None):
        '''Return the irreducible k-points corresponding to mesh

        Args:
            kgrid (int array): the kpoint grid
        '''
        try:
            assert np.shape(kgrid) == (3,)
            assert str(np.array(kgrid).dtype).startswith('int')
        except:
            raise SymmetryError("Invalid kgrid input: {}".format(kgrid))
        _mapping, _grid = spg.get_ir_reciprocal_mesh(kgrid, self._cell, is_shift=shift)
        
    
    def get_primitive(self):
        '''Return the primitive cell of the input cell.

        If primitive cell is not found, the original cell is returned

        Note:
            kwargs, except ``unit`` and ``coordSys`` are lost, when the primitive cell
            is found and returned.


        Returns:
            None, if the primitive cell is not found. True, if the original cell is already
            pritmitive, otherwise False.

            Cell or its subclass instance, depending on the ``latt`` at instantialization
        '''
        _flag = None
        _type = type(self._cell)
        _typeMap = self._cell.typeMapping
        _latt, _pos, _indice = spg.find_primitive(self._cell.get_spglib_input(), \
            symprec=self._symprec)
        _atoms = [_typeMap[v] for _i, v in enumerate(_indice)]
        _primCell = _type(_latt, _atoms, _pos, \
            unit=self._cell.unit, coordSys=self._cell.coordSys)
        if _primCell != None:
            # determine if the original cell is primitive
            if _primCell.vol == self._cell.vol:
                _flag = True
            else:
                _flag = False
            # the space group does not change when converting to primitive cell
            self._cell = _primCell
        return _flag, self._cell
