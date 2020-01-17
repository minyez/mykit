# coding = utf-8
'''class and functions related to crystal symmetry
'''

from collections import OrderedDict

import numpy as np
import spglib

from mykit.core.kmesh import (_check_valid_kpath_dict,
                              _check_valid_ksym_coord_pair, kpath_decoder)
from mykit.core.metadata._spk import _special_kpoints
from mykit.core.numeric import Prec


class SymmetryError(Exception):
    pass


def get_spacegroup(cell):
    '''return the space group number and symbol

    Args:
        cell (``Cell``)
    '''
    _spglib_check_cell_and_coordSys(cell)
    spg = spglib.get_spacegroup(cell.get_spglib_input(),
                                symprec=Prec._symprec).split()
    return int(spg[1][1:-1]), spg[0]


def get_sym_operations(cell):
    '''return symmetry operations.

    Args:
        cell (``Cell``)
    '''
    _spglib_check_cell_and_coordSys(cell)
    ops = spglib.get_symmetry(cell.get_spglib_input(), symprec=Prec._symprec)
    return ops


def get_ibzkpt(cell, kgrid, shift=None, frac=True):
    '''Return the irreducible k-points (IBZK) corresponding to mesh

    Args:
        cell (``Cell``)
        kgrid (int array): the kpoint grid
        frac (bool): if True, the fractional coordinate will be returned,
            otherwise the integer grid

    Returns:
        Two arrays, with n the number of IBZK
            - the coordinate of IBZK, shape (n,3) 
            - the un-normalized weight of each IBZK, shape (n,)
    '''
    _spglib_check_cell_and_coordSys(cell)
    _kgrid = np.array(kgrid)
    _mapping, _grid = spglib.get_ir_reciprocal_mesh(_kgrid,
                                                    cell.get_spglib_input(), is_shift=shift)

    _div = {True: _kgrid, False: np.array([1, 1, 1])}[frac]
    # All k-points and mapping to ir-grid points
    # ind_ibzk = list(set(_mapping))
    ind_ibzk = list(OrderedDict.fromkeys(_mapping, 0).keys())
    ibzk = []
    for i in ind_ibzk:
        ibzk.append(_grid[i, :]/_div)
    ibzk = np.array(ibzk, dtype=Prec._dtype)
    weight = OrderedDict.fromkeys(ind_ibzk, 0)
    for _i, v in enumerate(_mapping):
        weight[v] += 1
    return ibzk, np.array(list(weight.values()))


def get_inequiv(cell):
    '''Return the inequivalent atoms and there mapping
    '''
    _spglib_check_cell_and_coordSys(cell)
    data = spglib.get_symmetry_dataset(cell.get_spglib_input(), symprec=Prec._symprec)
    ea = data["equivalent_atoms"]
    typeMap = cell.typeMapping
    d = dict(((x, []) for x in set(ea)))
    for i, x in enumerate(ea):
        d[x].append(i)
    atomsEquiv = [typeMap[x] for x in d.keys()]
    indices = list(d.values())
    return atomsEquiv, indices


def is_primitive(cell):
    '''Check if the cell is primitive by its volume

    Args:
        cell (``Cell``)

    Returns:
        None, if the primitive cell is not found, often due to bad input
        True, if the input cell is primitve, False otherwise
    '''
    _spglib_check_cell_and_coordSys(cell)
    flag = None
    primCell = spglib.find_primitive(cell.get_spglib_input(),
                                     symprec=Prec._symprec)
    if primCell != None:
        vol = np.linalg.det(primCell[0])
        if abs(vol - cell.vol) < 0.1:
            flag = True
        else:
            flag = False
    return flag


def primitize(cell):
    '''Return the primitive cell of the input cell.

    If primitive cell is not found, the original cell is returned

    Args:
        cell (``Cell``)

    Returns:
    '''
    _spglib_check_cell_and_coordSys(cell)
    t = type(cell)
    typeMap = cell.typeMapping
    primCell = spglib.find_primitive(cell.get_spglib_input(),
                                     symprec=Prec._symprec)
    if primCell != None:
        latt, pos, indice = primCell
        atoms = [typeMap[v] for _i, v in enumerate(indice)]
        newCell = t(latt, atoms, pos, **cell.get_kwargs())
    else:
        newCell = cell
    return newCell


def standardize(cell, primitive=False):
    '''Return the standardized cell from the original cell

    Args:
        primitive (bool): if set True, the standard primitive 
            cell will be returned

    Returns:
        Cell or its subclass instance, depending on the ``cell`` 
        at instantialization
    '''
    assert isinstance(primitive, bool)
    _spglib_check_cell_and_coordSys(cell)
    
    t = type(cell)
    typeMap = cell.typeMapping
    stdCell = spglib.standardize_cell(cell.get_spglib_input(),
                                      to_primitive=primitive, symprec=Prec._symprec)
    if stdCell != None:
        latt, pos, indice = stdCell
        atoms = [typeMap[v] for _i, v in enumerate(indice)]
        newCell = t(latt, atoms, pos, **cell.get_kwargs())
    else:
        newCell = cell
    return newCell


# pylint: disable=bad-whitespace
class SpaceGroup:
    '''the space group symbols documented in the international table of crystallography A (ITA)

    TODO:
        check the symbols with spglib, especially for the underlines
    '''

    symbols = (
        'P1', 'P-1', 'P2', 'P21', 'C2',  # 1
        'Pm', 'Pc', 'Cm', 'Cc', 'P2m',
        'P21m', 'C2m', 'P2c', 'P21c', 'C2c',  # 11
        'P222', 'P2221', 'P21212', 'P212121', 'C2221',
        'C222', 'F222', 'I222', 'I212121', 'Pmm2',  # 21
        'Pmc21', 'Pcc2', 'Pma2', 'Pca21', 'Pnc2',
        'Pmn21', 'Pba2', 'Pna21', 'Pnn2', 'Cmm2',  # 31
        'Cmc21', 'Ccc2', 'Amm2', 'Aem2', 'Ama2',
        'Aea2', 'Fmm2', 'Fdd2', 'Imm2', 'Iba2',  # 41
        'Ima2', 'Pmmm', 'Pnnn', 'Pccm', 'Pban',
        'Pmma', 'Pnna', 'Pmna', 'Pcca', 'Pbam',  # 51
        'Pccn', 'Pbcm', 'Pnnm', 'Pmmn', 'Pbcn',
        'Pbca', 'Pnma', 'Cmcm', 'Cmce', 'Cmmm',  # 61
        'Cccm', 'Cmme', 'Ccce', 'Fmmm', 'Fddd',
        'Immm', 'Ibam', 'Ibca', 'Imma', 'P4',  # 71
        'P41', 'P42', 'P43', 'I4', 'I41',
        'P-4', 'I-4', 'P4m', 'P42m', 'P4n',  # 81
        'P42n', 'I4m', 'I41a', 'P422', 'P4212',
        'P4122', 'P41212', 'P4222', 'P42212', 'P4322',  # 91
        'P43212', 'I422', 'I4122', 'P4mm', 'P4bm',
        'P42cm', 'P42nm', 'P4cc', 'P4nc', 'P42mc',  # 101
        'P42bc', 'I4mm', 'I4cm', 'I41md', 'I41cd',
        'P-42m', 'P-42c', 'P-421m', 'P-421c', 'P-4m2',  # 111
        'P-4c2', 'P-4b2', 'P-4n2', 'I-4m2', 'I-4c2',
        'I-42m', 'I-42d', 'P4mmm', 'P4mcc', 'P4nbm',  # 121
        'P4nnc', 'P4mbm', 'P4mnc', 'P4nmm', 'P4ncc',
        'P42mmc', 'P42mcm', 'P42nbc', 'P42nnm', 'P42mbc',  # 131
        'P4_2/mnm', 'P42nmc', 'P42ncm', 'I4mmm', 'I4mcm',
        'I4_1/amd', 'I41acd', 'P3', 'P31', 'P32',  # 141
        'R3', 'P-3', 'R-3', 'P312', 'P321',
        'P3112', 'P3121', 'P3212', 'P3221', 'R32',  # 151
        'P3m1', 'P31m', 'P3c1', 'P31c', 'R3m',
        'R3c', 'P-31m', 'P-31c', 'P-3m1', 'P-3c1',  # 161
        'R-3m', 'R-3c', 'P6_3mc', 'P61', 'P65',
        'P62', 'P64', 'P63', 'P-6', 'P6m',  # 171
        'P63m', 'P622', 'P6122', 'P6522', 'P6222',
        'P6422', 'P6322', 'P6mm', 'P6cc', 'P63cm',  # 181
        'P6_3mc', 'P-6m2', 'P-6c2', 'P-62m', 'P-62c',
        'P6mmm', 'P6mcc', 'P63mcm', 'P63mmc', 'P23',  # 191
        'F23', 'I23', 'P213', 'I2_13', 'Pm-3',
        'Pn-3', 'Fm-3', 'Fd-3', 'Im-3', 'Pa-3',  # 201
        'Ia-3', 'P432', 'P4232', 'F432', 'F4132',
        'I432', 'P4332', 'P4132', 'I4132', 'P-43m',  # 211
        'F-43m', 'I-43m', 'P-43n', 'F-43c', 'I-43d',
        'Pm-3m', 'Pn-3n', 'Pm-3n', 'Pn-3m', 'Fm-3m',  # 221
        'Fm-3c', 'Fd-3m', 'Fd-3c', 'Im-3m', 'Ia-3d',
    )

    @staticmethod
    def k_trans_mat_from_prim_to_conv(iden):
        '''Return transformation matrix to convert k-point coordinate in 
        reciprocal lattice vectors of primitive cell to that in vectors of conventional cell

        The transormation matrix is retrived from Bilbao server.
        The k-point coordinate in primitive cell, (u,v,w), is converted to 
        conventional cell (p,q,n) by the matrix A as

        (p,q,n)^T = A x (u,v,w)^T

        i.e. both coordinate are treated as a column vector

        Note:
            For space group 38,39,40,41, tranformation matrix to conventional basis C2mm is used

        Args:
            iden (int): the id of space group of the cell, 1~230

        Returns:
            array: (3,3)
        '''
        _check_valid_spg_id(iden)
        # primitive and conventional cells coincide.
        identity = np.array(
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])

        _dict = {
            (5, 8, 9, 12, 15,):
            np.array([[1.0, -1.0, 0.0],
                      [1.0, 1.0, 0.0],
                      [0.0, 0.0, 1.0]]),  # u-v, u+v, w
            (20, 21, 35, 36, 37, 38, 39, 40, 41, 63, 64, 65, 66, 67, 68,):
            np.array([[1.0, 1.0, 0.0],
                      [-1.0, 1.0, 0.0],
                      [0.0, 0.0, 1.0]]),  # u+v, -u+v, w
            (22, 42, 43, 69, 70, 196, 202, 203, 209, 210, 216, 219, 225,\
             226, 227, 228):
            np.array([[-1.0, 1.0, 1.0],
                      [1.0, -1.0, 1.0],
                      [1.0, 1.0, -1.0]]),  # -u+v+w, u-v+w, u+v-w
            (23, 24, 44, 45, 46, 71, 72, 73, 74, 79, 80, 82, 87, 88, 97, 98,\
             107, 108, 109, 110, 119, 120, 121, 122, 139, 140, 141, 142,\
             197, 204, 206, 211, 214, 217, 220, 229, 230):
            np.array([[0.0, 1.0, 1.0],
                      [1.0, 0.0, 1.0],
                      [1.0, 1.0, 0.0]]),  # v+w, u+w, u+v
            (146, 148, 155, 160, 161, 166, 167):
            np.array([[1.0, -1.0, 0.0],
                      [0.0, 1.0, -1.0],
                      [1.0, 1.0, 1.0]]),  # u-v, v-w, u+v+w
        }
        for k, v in _dict.items():
            if iden in k:
                return v
        return identity

    @classmethod
    def get_spg_id(cls, symbol):
        '''Get the id in ITA from the symbol of space group of symbol

        Args:
            symbol (str): the symbol of space group
        '''
        if isinstance(symbol, int):
            raise SymmetryError("int received. symbol should be string.")
        try:
            _id = 1 + cls.symbols.index(symbol)
        except ValueError:
            raise SymmetryError(
                "Space group symbol is not found in ITA: {}".format(symbol))
        return _id

    @classmethod
    def get_spg_symbol(cls, iden):
        '''Get the symbol of space group with index ``id`` in ITA

        Args:
            iden (int): the id of space group, 1~230
        '''
        _check_valid_spg_id(iden)
        return cls.symbols[iden-1]


class SpecialKpoints(Prec):
    '''class for special kpoints of space groups

    Args:
        iden (int): the id of space group
        alen (array): (3,) array of lattice constant a,b,c of conventional cell
        isPrimitive: if set True, coordinates in primitive cell
        will be used.
        custom_symbols (dict): dictionary of custom kpoints, each item
        being a symbol-coordinate pair. The coordinate must be in
        the system of primitive cell.
        Note that this will overwrite the default definition.
    '''

    def __init__(self, iden, alen, isPrimitive, custom_symbols=None):
        _check_valid_spg_id(iden)
        try:
            assert np.shape(alen) == (3,)
        except AssertionError:
            raise SymmetryError(
                "Invalid lattice constant shape: {}".format(alen))

        self.spgId = iden
        self._spDict = _special_kpoints[iden]
        try:
            iset = self._spDict["cond"](*alen)
        except ValueError as _err:
            raise SymmetryError(str(_err))
        self._sp = self._spDict["spPrim"][iset]
        self._isPrim = isPrimitive
        self._transMat = SpaceGroup.k_trans_mat_from_prim_to_conv(self.spgId)
        self._custom = _check_valid_custom_ksym_dict(custom_symbols)
        try:
            self._kpaths = self._spDict["kpath"][iset]
        except KeyError:
            self._kpaths = []

    def __getitem__(self, ksym):
        '''Get the coordinate of kpoint symbol.

        Custom symbols will be searched first.
        '''
        coord = self._custom.get(ksym, None)
        if coord is None:
            coord = self._sp.get(ksym, None)
        if coord is None:
            raise IndexError("kpoint symbol not defined: {}".format(ksym))
        if not self._isPrim:
            coord = np.dot(coord, np.transpose(self._transMat))
        return list(coord)

    def __setitem__(self, ksym, coord):
        _check_valid_ksym_coord_pair(ksym, coord)
        self._custom.update({ksym: coord})

    @property
    def spkSym(self):
        '''List. Symbols of all available special kpoints
        '''
        return list(self.spkCoord.keys())

    @property
    def spkCoord(self):
        '''List. Coordinates in primitive cell of all available special kpoints
        '''
        _ret = {}
        _ret.update(self._sp)
        _ret.update(self._custom)
        return _ret

    def check_kpaths_predef(self, ipath=None):
        '''Return the predefined kpath string

        Args:
            ipath (int): the index of predefined kpath string to check.
            None to return all kpath strings in a list.
            If ipath is not found, an empty string will be returned.

        Returns:
            list if ipath is None, otherwise string
        '''
        kpath = self._kpaths
        if ipath != None:
            try:
                assert isinstance(ipath, int)
            except AssertionError:
                raise ValueError("ipath should be int")
            try:
                kpath = self._kpaths[ipath]
            except IndexError:
                kpath = ''
        return kpath

    def convert_kpath(self, kpathStr):
        '''Convert a kpath string to coordinates

        Args:
            kpathStr (str): the kpath to translate

        Returns:
            dict with two keys, namely "symbols" and "coordinates"
        '''
        try:
            assert isinstance(kpathStr, str)
        except AssertionError:
            raise ValueError(
                "input kpath should be a str, not {}".format(type(kpathStr)))
        _ret = {}
        decodedSyms = kpath_decoder(kpathStr)
        coords = []
        for ksym in decodedSyms:
            if ksym in self._custom:
                coords.append(self._custom[ksym])
            elif ksym in self.spkCoord:
                coords.append(self.spkCoord[ksym])
            else:
                raise SymmetryError(
                    "kpoint symbol not defined: {}".format(ksym))
        coords = np.array(coords, dtype=self._dtype)
        if not self._isPrim:
            coords = np.dot(coords, np.transpose(self._transMat))
        _ret["symbols"] = decodedSyms
        _ret["coordinates"] = coords
        # self-check
        _check_valid_kpath_dict(_ret)
        return _ret

    def convert_kpaths_predef(self, ipath=None):
        '''Get the coordinates of special kpoints along predefined kpath

        Args:
            ipath (int): the index of predefined kpath string in check_kpaths_predef
        '''
        kpathPredef = self.check_kpaths_predef(ipath)
        if not kpathPredef in ([], ''):
            if isinstance(kpathPredef, str):
                return self.convert_kpath(kpathPredef)
            if isinstance(kpathPredef, list):
                return list(map(self.convert_kpath, kpathPredef))
        return None

    # @classmethod
    # def from_symmetry(cls, sym, custom_symbols=None):
    #     '''Create special kpoints instance from a Symmetry instance

    #     Note:
    #         It is safer to use sym with its cell standardized,
    #         otherwise the condition check for special kpoints set
    #         may fail.

    #     Args:
    #         sym: Symmetry instance
    #     '''
    #     try:
    #         assert isinstance(sym, Symmetry)
    #     except:
    #         raise SymmetryError("The input should be Symmetry instance.")
    #     # when primitive cell is meeted, need to extract the
    #     # lattice constants of conventional or standard cell
    #     isPrim = sym.isPrimitive
    #     if isPrim is None:
    #         raise SymmetryError(
    #             "Fail to determine whether the cell is primitive. Check the cell")
    #     elif isPrim:
    #         flag, stdCell = sym.get_standard()
    #         # I think the exception is redudant, since if the cell cannot
    #         # be standardized, isPrim flag is None and should be captured
    #         # above. But anyway it is put here for safety.
    #         if not flag:
    #             raise SymmetryError(
    #                 "Fail to standardize the cell to obtain lattice constants")
    #         alen = stdCell.alen
    #     else:
    #         alen = sym.alen
    #     return cls(sym.spgId, alen, isPrim, custom_symbols=custom_symbols)

    @classmethod
    def from_cell(cls, cell, custom_symbols=None):
        '''Create special kpoints instance from a Cell instance.

        Note:
            It is safer to use a standardized cell, otherwise the 
            condition check for special kpoints set may fail.

        Args:
            cell (``Cell``)
            custom_symbols (dict): custom symbol-coordinate pair
        '''
        isPrim = is_primitive(cell)
        spgId, _ = get_spacegroup(cell)
        if isPrim is None:
            raise SymmetryError(
                "Fail to determine whether the cell is primitive. Check the cell")
        elif isPrim:
            stdCell = standardize(cell, primitive=False)
            alen = stdCell.alen
        else:
            alen = cell.alen
        return cls(spgId, alen, isPrim, custom_symbols=custom_symbols)

    @classmethod
    def get_kpath_from_cell(cls, pathStr, cell, custom_symbols=None):
        '''Return coordinates of all predefined kpaths for the space group of the cell

        Args:
            pathStr (str)
            cell (``Cell``)
            custom_symbols (dict): custom symbol-coordinate pair
        '''
        spk = cls.from_cell(cell, custom_symbols=custom_symbols)
        return spk.convert_kpath(pathStr)

    @classmethod
    def get_kpaths_predef_from_cell(cls, cell, ipath=None, custom_symbols=None):
        '''Return coordinates of all predefined kpaths for the space group of the cell

        Args:
            cell: instance of Cell or its subclass
            ipath (int): the index of predefined path.
                If not specified, all paths will be returned
            custom_symbols (dict): custom symbol-coordinate pair
        '''
        _spglib_check_cell_and_coordSys(cell)
        spk = cls.from_cell(cell, custom_symbols=custom_symbols)
        return spk.convert_kpaths_predef(ipath=ipath)


def _spglib_check_cell_and_coordSys(cellIn):
    '''Raise if the input is not an instance of Cell or its subclasses,
    or its coordinate system is not direct (i.e. fractional), required by spglib

    Args:   
        cellIn (any type): the input to check whether it is a cell instance
    '''
    from mykit.core.cell import Cell
    try:
        assert isinstance(cellIn, Cell)
    except AssertionError:
        raise SymmetryError(
            "The input should be instance of Cell or its subclass.")
    try:
        assert cellIn.coordSys == "D"
    except AssertionError:
        raise SymmetryError(
            "The coordinate system should be direct. Cartisian found.")


def _check_valid_spg_id(iden):
    '''Raise if iden is not a valid spacegroup id, i.e. 1~230

    Args:
        iden (int): the space group id to check
    '''
    if isinstance(iden, str):
        raise SymmetryError("string received. id should be int.")
    if not isinstance(iden, int):
        raise SymmetryError("id should be int.")
    try:
        assert iden in range(1, 231)
    except AssertionError:
        raise SymmetryError("Invalid space group id (1~230): {}".format(iden))


def _check_valid_custom_ksym_dict(custom_symbols):
    '''check if a custom symbol dictionary is valid
    '''
    if custom_symbols != None:
        try:
            assert isinstance(custom_symbols, dict)
        except:
            raise TypeError(
                "custom_symbols must be dict, received {}".format(type(custom_symbols)))
        for k, v in custom_symbols.items():
            _check_valid_ksym_coord_pair(k, v)
        return custom_symbols
    else:
        return {}
