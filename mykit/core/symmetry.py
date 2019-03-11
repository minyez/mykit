# coding = utf-8
'''class and functions related to crystal symmetry
'''

import numpy as np
import spglib as spg

from mykit.core.lattice import lattice
from mykit.core.numeric import prec


class SymmetryError(Exception):
    pass


class symmetry(prec):
    '''the class for symmetry information of crystal, powered by spglib

    Args:
        latt (lattice or its subclasses)
    '''
    
    def __init__(self, latt):

        try:
            assert isinstance(latt, lattice)
        except AssertionError:
            raise SymmetryError("The input should be instance of lattice or its subclass.")
        try:
            assert latt.unit == "D"
        except AssertionError:
            raise SymmetryError("The coordiante system should be direct. Cartisian found.")
            
        # convert to direct coordinate system
        self._lattice = latt
        self.spaceGroup = spg.get_spacegroup(self._lattice.get_spglib_input(), \
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
        _mapping, _grid = spg.get_ir_reciprocal_mesh(kgrid, self._lattice, is_shift=shift)
        
    
    def get_primitive(self):
        '''Return the primitive cell of the input lattice.

        If primitive cell is not found, the original lattice is returned

        Note:
            kwargs, except ``unit`` and ``coordSys`` are lost, when the primitive cell
            is found and returned.


        Returns:
            None, if the primitive cell is not found. True, if the original cell is already
            pritmitive, otherwise False.

            lattice or its subclass instance, depending on the ``latt`` at instantialization
        '''
        _flag = None
        _type = type(self._lattice)
        _typeMap = self._lattice.typeMapping
        _cell, _pos, _indice = spg.find_primitive(self._lattice.get_spglib_input(), \
            symprec=self._symprec)
        _atoms = [_typeMap[v] for _i, v in enumerate(_indice)]
        _primLatt = _type(_cell, _atoms, _pos, \
            unit=self._lattice.unit, coordSys=self._lattice.coordSys)
        if _primLatt != None:
            # determine if the original lattice is primitive
            if _primLatt.vol == self._lattice.vol:
                _flag = True
            else:
                _flag = False
            # the space group does not change when converting to primitive cell
            self._lattice = _primLatt
        return _flag, self._lattice
