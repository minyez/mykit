# coding = utf-8
import numpy as np

from mykit.core.log import Verbose
from mykit.core.numeric import Prec
from mykit.core.unit import EnergyUnit
from mykit.core.utils import get_str_indices_by_iden

KEYS_DOS_PROJ = ("atoms", "projs", "pDos")
'''Tuple. Required keys for projected DOS input

"atoms": a list of strings

"projs": a list of strings

"pDos": array-like, floats as elements. 
The shape should be (nedos, nspins, natoms, nprojs)
'''


class DosError(Exception):
    pass


class Dos(Prec, Verbose, EnergyUnit):
    '''Base class for analyzing density of states (DOS) data

    The energy grids of DOS, total DOS and Fermi level are required.
    The energy grids should be 1d array-like object, with ``nedos`` as length.
    The shape of total DOS should be (nedos, nspins)

    Args:
        edos (1d-array-like)
        totalDOS (2d-array-like)
        efermi (float)

    Optional args:
        unit ('ev','ry','au'): the unit of the energy grid points, in lower case.
        projected (dict): projected DOS.
            It should have three keys, "atoms", "projs" and "pDos"
    '''

    def __init__(self, edos, totalDos, efermi, unit='ev', projected=None):
        # check shape consistency
        try:
            shapeE = np.shape(edos)
            shapeDos = np.shape(totalDos)
            assert len(shapeE) == 1
            assert len(shapeDos) == 2
            assert shapeDos[0] == shapeE[0]
        except:
            raise DosError(
                'Inconsistent shape: edos {}, DOS {}'.format(shapeE, shapeDos))
        EnergyUnit.__init__(self, eunit=unit)
        self._edos = np.array(edos, dtype=self._dtype)
        self._efermi = efermi
        self._nedos, self._nspins = shapeDos
        self._dos = np.array(totalDos, dtype=self._dtype)
        self._projParsed = False
        self._parse_proj(projected)

    @property
    def unit(self):
        return self._eunit

    @unit.setter
    def unit(self, newu):
        coef = self._get_eunit_conversion(newu)
        if coef != 1:
            self._edos *= coef
            self._efermi *= coef
            self._eunit = newu.lower()

    @property
    def edos(self):
        '''float array. Energy grid points of DOS'''
        return self._edos

    @property
    def nedos(self):
        '''int. Number of energy grid points'''
        return self._nedos

    @property
    def nspins(self):
        '''int. Number of spin channels'''
        return self._nspins

    @property
    def dos(self):
        '''float array. Total density of states

        Shape (nedos, nspins)
        '''
        return self._dos

    @property
    def efermi(self):
        '''float. Fermi level'''
        return self._efermi

    def _parse_proj(self, projected):
        '''parse the projected DOS information
        '''
        if not projected is None:
            try:
                self._atoms, self._projs, pDos = \
                    self._check_project_consistency(projected)
                self._pDos = np.array(pDos, dtype=self._dtype)
                self._projParsed = True
            except ValueError:
                self.print_warn("Bad projection input. Skip.")

    def _check_project_consistency(self, projected):
        try:
            assert isinstance(projected, dict)
        except AssertionError:
            return ()
        try:
            for key in KEYS_DOS_PROJ:
                assert key in projected
        except AssertionError:
            return ()
        # print(projected)
        atoms = projected["atoms"]
        projs = projected["projs"]
        pDos = projected["pDos"]
        try:
            natoms = len(atoms)
            nprojs = len(projs)
            self.print_log("Shapes: ", np.shape(pDos),
                           (self._nedos, self._nspins,
                            natoms, nprojs),
                           level=3)
            assert np.shape(pDos) == \
                (self._nedos, self._nspins, natoms, nprojs)
        except (AssertionError, TypeError):
            return ()
        return atoms, projs, pDos

    @property
    def hasProjection(self):
        '''Bool, whether the ``Dos`` object has projected DOS'''
        return self._projParsed

    @property
    def atoms(self):
        '''list of string. The symbols of atoms
        on which the DOS is projected

        None if no projection was parsed
        '''
        if self.hasProjection:
            return self._atoms
        return None

    @property
    def natoms(self):
        '''int. Number of atoms

        0 if no projection was parsed
        '''
        try:
            return len(self.atoms)
        except TypeError:
            return 0

    @property
    def projs(self):
        '''list of string. The symbols of orbital projectors
        on which the DOS is projected

        None if no projection was parsed
        '''
        if self.hasProjection:
            return self._projs
        return None

    @property
    def nprojs(self):
        '''int. Number of projectors

        0 if no projection was parsed
        '''
        try:
            return len(self.projs)
        except TypeError:
            return 0

    @property
    def pDos(self):
        '''float array. Projected DOS

        shape (nedos, nspins, natoms, nprojs)
        '''
        if self.hasProjection:
            return self._pDos
        return None

    def sum_atom_proj_comp(self, atom=None, proj=None, fail_one=True):
        '''Sum the pDOS for projectors `proj` on atoms `atom`

        If the instance does not have projection information, 1 will be returned.

        Args:
            atom (int, str, Iterable)
            proj (int, str, Iterable)
            fail_one (bool): control the return when no projection is available. 
                if set True, return np.ones with correct shape, otherwise np.zeros

        Returns:
            (nedos, nspins)
        '''
        if not self.hasProjection:
            func = {True: np.ones, False: np.zeros}
            try:
                return func[fail_one]((self.nedos, self.nspins))
            except KeyError:
                raise TypeError("fail_one should be bool type.")
        if atom is None:
            atInd = list(range(self.natoms))
        else:
            atInd = self._get_atom_indices(atom)
        if proj is None:
            prInd = list(range(self.nprojs))
        else:
            prInd = self._get_proj_indices(proj)
        coeff = np.zeros((self.nedos, self.nspins))
        for a in atInd:
            for p in prInd:
                coeff += self.pDos[:, :, a, p]
        return coeff

    def _get_atom_indices(self, atom):
        if self.hasProjection:
            return get_str_indices_by_iden(self._atoms, atom)
        return []

    def _get_proj_indices(self, proj):
        if self.hasProjection:
            return get_str_indices_by_iden(self._projs, proj)
        return []
