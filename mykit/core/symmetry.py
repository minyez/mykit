# coding = utf-8
'''class and functions related to crystal symmetry
'''

import json
import os
from collections import OrderedDict

import numpy as np
import spglib

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
        _spg = spglib.get_spacegroup(self._cell.get_spglib_input(), \
            symprec=self._symprec).split()
        self._spgSym = _spg[0]
        self._spgId = int(_spg[1][1:-1])

    @property
    def operations(self):
        '''Symmetry operations. See spglib get_symmetry docstring
        '''
        _ops = spglib.get_symmetry(self._cell, symprec=self._symprec)
        return [(r, t) for r, t in zip(_ops["rotations"], _ops["translations"])]

    @property
    def spgSym(self):
        return self._spgSym

    @property
    def spgId(self):
        return self._spgId

    def ibzkpt(self, kgrid, shift=None):
        '''Return the irreducible k-points (IBZK) corresponding to mesh

        Args:
            kgrid (int array): the kpoint grid
        
        Returns:
            Two arrays, with n the number of IBZK
            (n,3) the fractional coordinate of IBZK in reciprocal lattice vectors
            (n,) the un-normalized weight of each IBZK
        '''
        # * leave the type check to spglib
        # try:
        #     assert np.shape(kgrid) == (3,)
        #     assert str(np.array(kgrid).dtype).startswith('int')
        # except AssertionError:
        #     raise SymmetryError("Invalid kgrid input: {}".format(kgrid))
        _kgrid = np.array(kgrid)
        _mapping, _grid = spglib.get_ir_reciprocal_mesh(_kgrid, \
            self._cell.get_spglib_input(), is_shift=shift)
   
        # All k-points and mapping to ir-grid points
        ind_ibzk = list(set(_mapping))
        ibzk = []
        for i in ind_ibzk:
            ibzk.append(_grid[i, :]/_kgrid)
        ibzk = np.array(ibzk, dtype=self._dtype)
        weight = OrderedDict.fromkeys(ind_ibzk, 0)
        for _i, v in enumerate(_mapping):
            weight[v] += 1
        return ibzk, np.array(list(weight.values()))
        
    def get_primitive(self):
        '''Return the primitive cell of the input cell.

        If primitive cell is not found, the original cell is returned

        Note:
            kwargs, except ``unit`` and ``coordSys`` are lost, when the primitive cell
            is found and returned.


        Returns:
            None, if the primitive cell is not found. True, if the original cell is already
            primititive, otherwise False.

            Cell or its subclass instance, depending on the ``cell`` at instantialization
        '''
        _flag = None
        _type = type(self._cell)
        _typeMap = self._cell.typeMapping
        _primCell = spglib.find_primitive(self._cell.get_spglib_input(), \
            symprec=self._symprec)
        if _primCell != None:
            _latt, _pos, _indice = _primCell
            _atoms = [_typeMap[v] for _i, v in enumerate(_indice)]
            _newCell = _type(_latt, _atoms, _pos, \
                unit=self._cell.unit, coordSys=self._cell.coordSys)
            # determine if the original cell is primitive
            if _newCell.vol == self._cell.vol:
                _flag = True
            else:
                _flag = False
        else:
            _newCell = self._cell
        return _flag, _newCell

    def get_standard(self, primitive=False):
        '''Return the standardized cell from the original cell

        Args:
            primitive (bool): if set True, the standard primitive cell will be returned
        
        Returns:
            Cell or its subclass instance, depending on the ``cell`` at instantialization
        '''
        assert isinstance(primitive, bool)
        _type = type(self._cell)
        _typeMap = self._cell.typeMapping
        _stdCell = spglib.standardize_cell(self._cell.get_spglib_input(), \
            to_primitive=primitive, symprec=self._symprec)
        if _stdCell != None:
            _latt, _pos, _indice = _stdCell
            _atoms = [_typeMap[v] for _i, v in enumerate(_indice)]
            _newCell = _type(_latt, _atoms, _pos, \
                unit=self._cell.unit, coordSys=self._cell.coordSys)
        else:
            _newCell = self._cell
        return _newCell


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
    
    @staticmethod
    def k_trans_mat_from_prim_to_conv(id):
        '''Convert k-point coordinate in reciprocal lattice of primitive cell to that of conventional cell

        The transormation matrix is retrived from Bilbao server.
        The k-point coordinate in primitive cell, (u,v,w), is converted to conventional cell (p,q,n) by
        the matrix A as

        (p,q,n)^T = A x (u,v,w)^T

        i.e. both coordinate are treated as a column vector

        Note:
            For space group 38,39,40,41, tranformation matrix to conventional basis C2mm is used

        Args:
            id (int): the id of space group, 1~230

        Returns:
            array: (3,3)
        '''
        try:
            assert isinstance(id, int)
        except AssertionError:
            raise SymmetryError("Invalid space group id (1~230): {}".format(id))
        identity = np.array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]])
        _dict = {
            (5,8,9,12,15):
                   np.array([[ 1.0,-1.0, 0.0],
                             [ 1.0, 1.0, 0.0],
                             [ 0.0, 0.0, 1.0]]), # u-v, u+v, w
            (20,21,35,36,37,38,39,40,41):
                   np.array([[ 1.0, 1.0, 0.0],
                             [-1.0, 1.0, 0.0],
                             [ 0.0, 0.0, 1.0]]), # u+v, -u+v, w
            (22,42,43):
                   np.array([[-1.0, 1.0, 1.0],
                             [ 1.0,-1.0, 1.0],
                             [ 1.0, 1.0,-1.0]]), # -u+v+w, u-v+w, u+v-w
            (23,24,44,):
                   np.array([[ 0.0, 1.0, 1.0],
                             [ 1.0, 0.0, 1.0],
                             [ 1.0, 1.0, 0.0]]), # v+w, u+w, u+v
        }
        for k, v in _dict.items():
            if id in k:
                return v
        return identity



    @classmethod
    def get_spg_index(cls, symbol):
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
    def get_spg_symbol(cls, id):
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


class special_kpoints:
    '''class for special kpoints of space groups
    '''
    _meta = os.path.join(os.path.dirname(__file__), 'metadata', "special_kpoints.json")
    with open(_meta, 'r') as h:
        _j = json.load(h)

    def __init__(self, id):
        pass
