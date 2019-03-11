# coding = utf-8
'''class and functions related to crystal symmetry
'''

import numpy as np
import spglib as spg

from mykit.core.cell import Cell
from mykit.core.numeric import prec


class SymmetryError(Exception):
    pass


class Symmetry(prec):
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
            assert cell.coordSys == "D"
        except AssertionError:
            raise SymmetryError("The coordiante system should be direct. Cartisian found.")
            
        # convert to direct coordinate system
        self._cell = cell
        _spaceg = spg.get_spacegroup(self._cell.get_spglib_input(), \
            symprec=self._symprec).split()
        self._spgSym = _spaceg[0]
        self._spgId = int(_spaceg[1][1:-1])
    
    @property
    def spgSym(self):
        return self._spgSym

    @property
    def spgId(self):
        return self._spgId

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


# pylint: disable=bad-whitespace
class space_group:
    '''the space group symbols documented in the international table of crystallography (ITC)
    '''

    symbols = (
            'P1'    , 'P-1'   , 'P2'    , 'P21'    , 'C2'    ,
            'Pm'    , 'Pc'    , 'Cm'    , 'Cc'     , 'P2m'   ,
            'P21m'  , 'C2m'   , 'P2c'   , 'P21c'   , 'C2c'   ,
            'P222'  , 'P2221' , 'P21212', 'P212121', 'C2221' ,
            'C222'  , 'F222'  , 'I222'  , 'I212121', 'Pmm2'  ,
            'Pmc21' , 'Pcc2'  , 'Pma2'  , 'Pca21'  , 'Pnc2'  ,
            'Pmn21' , 'Pba2'  , 'Pna21' , 'Pnn2'   , 'Cmm2'  ,
            'Cmc21' , 'Ccc2'  , 'Amm2'  , 'Aem2'   , 'Ama2'  ,
            'Aea2'  , 'Fmm2'  , 'Fdd2'  , 'Imm2'   , 'Iba2'  ,
            'Ima2'  , 'Pmmm'  , 'Pnnn'  , 'Pccm'   , 'Pban'  ,
            'Pmma'  , 'Pnna'  , 'Pmna'  , 'Pcca'   , 'Pbam'  ,
            'Pccn'  , 'Pbcm'  , 'Pnnm'  , 'Pmmn'   , 'Pbcn'  ,
            'Pbca'  , 'Pnma'  , 'Cmcm'  , 'Cmce'   , 'Cmmm'  ,
            'Cccm'  , 'Cmme'  , 'Ccce'  , 'Fmmm'   , 'Fddd'  ,
            'Immm'  , 'Ibam'  , 'Ibca'  , 'Imma'   , 'P4'    ,
            'P41'   , 'P42'   , 'P43'   , 'I4'     , 'I41'   ,
            'P-4'   , 'I-4'   , 'P4m'   , 'P42m'   , 'P4n'   ,
            'P42n'  , 'I4m'   , 'I41a'  , 'P422'   , 'P4212' ,
            'P4122' , 'P41212', 'P4222' , 'P42212' , 'P4322' ,
            'P43212', 'I422'  , 'I4122' , 'P4mm'   , 'P4bm'  ,
            'P42cm' , 'P42nm' , 'P4cc'  , 'P4nc'   , 'P42mc' ,
            'P42bc' , 'I4mm'  , 'I4cm'  , 'I41md'  , 'I41cd' ,
            'P-42m' , 'P-42c' , 'P-421m', 'P-421c' , 'P-4m2' ,
            'P-4c2' , 'P-4b2' , 'P-4n2' , 'I-4m2'  , 'I-4c2' ,
            'I-42m' , 'I-42d' , 'P4mmm' , 'P4mcc'  , 'P4nbm' ,
            'P4nnc' , 'P4mbm' , 'P4mnc' , 'P4nmm'  , 'P4ncc' ,
            'P42mmc', 'P42mcm', 'P42nbc', 'P42nnm' , 'P42mbc',
            'P42mnm', 'P42nmc', 'P42ncm', 'I4mmm'  , 'I4mcm' ,
            'I41amd', 'I41acd', 'P3'    , 'P31'    , 'P32'   ,
            'R3'    , 'P-3'   , 'R-3'   , 'P312'   , 'P321'  ,
            'P3112' , 'P3121' , 'P3212' , 'P3221'  , 'R32'   ,
            'P3m1'  , 'P31m'  , 'P3c1'  , 'P31c'   , 'R3m'   ,
            'R3c'   , 'P-31m' , 'P-31c' , 'P-3m1'  , 'P-3c1' ,
            'R-3m'  , 'R-3c'  , 'P6'    , 'P61'    , 'P65'   ,
            'P62'   , 'P64'   , 'P63'   , 'P-6'    , 'P6m'   ,
            'P63m'  , 'P622'  , 'P6122' , 'P6522'  , 'P6222' ,
            'P6422' , 'P6322' , 'P6mm'  , 'P6cc'   , 'P63cm' ,
            'P63mc' , 'P-6m2' , 'P-6c2' , 'P-62m'  , 'P-62c' ,
            'P6mmm' , 'P6mcc' , 'P63mcm', 'P63mmc' , 'P23'   ,
            'F23'   , 'I23'   , 'P213'  , 'I213'   , 'Pm-3'  ,
            'Pn-3'  , 'Fm-3'  , 'Fd-3'  , 'Im-3'   , 'Pa-3'  ,
            'Ia-3'  , 'P432'  , 'P4232' , 'F432'   , 'F4132' ,
            'I432'  , 'P4332' , 'P4132' , 'I4132'  , 'P-43m' ,
            'F-43m' , 'I-43m' , 'P-43n' , 'F-43c'  , 'I-43d' ,
            'Pm-3m' , 'Pn-3n' , 'Pm-3n' , 'Pn-3m'  , 'Fm-3m' ,
            'Fm-3c' , 'Fd-3m' , 'Fd-3c' , 'Im-3m'  , 'Ia-3d' ,
            )
    try:
        assert len(symbols) == 230
    except AssertionError:
        raise SymmetryError("Bad space group table. Contact developer.")


    @classmethod
    def get_index(cls, symbol):
        '''Get the index in ITC from the symbol of space group of symbol

        Args:
            symbol (str): the symbol of space group
        '''
        if isinstance(symbol, int):
            raise SymmetryError("int received. symbol should be string.")
        try:
            _id = 1 + cls.symbols.index(symbol)
        except ValueError:
            raise SymmetryError("Space group symbol is not found in ITC: {}".format(symbol))
        return _id
    
    @classmethod
    def get_symbol(cls, id):
        '''Get the symbol of space group with index ``id`` in ITC

        Args:
            id (int): the id of space group, 1~230
        '''
        if isinstance(id, str):
            raise SymmetryError("string received. id should be int.")
        try:
            assert isinstance(id, int)
            assert 0 < id < 231
        except AssertionError:
            raise SymmetryError("Invalid space group id (1~230): {}".format(id))
        return cls.symbols[id-1]
