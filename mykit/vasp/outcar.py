# coding = utf-8
'''Module that defines classes for VASP output files
'''
import os
import re
import subprocess as sp

import numpy as np

from mykit.core.bandstructure import BandStructure
from mykit.core.cell import atoms_from_sym_nat
from mykit.core.dos import Dos
from mykit.core.log import Verbose
from mykit.core.numeric import Prec
from mykit.core.utils import conv_string, trim_before
from mykit.vasp.incar import Incar
from mykit.vasp.poscar import Poscar


class OutcarError(Exception):
    pass


class Outcar(Verbose, Prec):

    def __init__(self, pathOut='OUTCAR'):

        # read in OUTCAR
        with open(pathOut, 'r') as f:
            self._outlines = f.readlines()

        self._efermi = None
        self._poscarComm = 'POSCAR Exported from OUTCAR'
        self._check_finished()
        self._read_all()
        # Build initial POSCAR
        self._initPoscar = Poscar(self._initLatt, self.atoms, self._initPos,
                                  unit="ang", coordSys="D", comment=self._poscarComm)

    @property
    def initPoscar(self):
        return self._initPoscar
    
    @property
    def efermi(self):
        '''float. Fermi energy for final structure'''
        if self._efermi is None:
            st, ed = self._get_ionic_data_region()
            while st <= ed:
                l = self._outlines[st].strip()
                if l.startswith('E-fermi'):
                    self._efermi = conv_string(l, float, 2)
                    break
                st += 1
        return self._efermi

    def _check_finished(self):
        l = self._outlines[-1].strip()
        if l.startswith('Voluntary context switches'):
            self._finished = True
        else:
            self._finished = False
        if not self._finished:
            self.print_warn(
                "Task has not finished. Try to load anyway...", level=0)

    def _read_all(self):
        '''initialize
        '''
        nl = len(self._outlines)
        i = 0
        reader = {
            r'INCAR:': self._read_atomtypes,
            r'Subroutine IBZKPT': self._read_weight,
            r'direct lattice vectors': self._read_init_latt,
            r'k-points in reciprocal lattice and weights': self._read_kpts,
            r'position of ions in fractional': self._read_init_pos,
            r'Dimension of arrays': self._read_dim_params,
            r'SYSTEM =': self._read_incar_params,
            r'k-point  1 :': self._read_npw,
            r'First call to EWALD': self._init_iterations,
        }
        while i < nl:
            l = self._outlines[i].strip().strip('-=')
            if l != '':
                for k, v in reader.items():
                    if re.match(k, l):
                        i = v(i)
            i += 1

    def _read_atomtypes(self, linenum):
        '''Reading the atom types from the POTCAR names

        Args:
            linenum (int): the index of starting line, 'INCAR:'

        Returns:
            int
        '''
        self._atomTypes = []
        first = self._outlines[linenum+1]
        i = 0
        while True:
            i += 1
            same = self._outlines[linenum + 1 + i] == first
            isPotcar = self._outlines[linenum +
                                      2 + i].strip().startswith('VRHFIN')
            if same and isPotcar:
                break
        ntypes = i
        for i in range(ntypes):
            symbol = self._outlines[linenum + 1 + i].split()[-2]
            self._atomTypes.append(symbol)
        return linenum + ntypes

    @property
    def atomTypes(self):
        return self._atomTypes

    def _read_init_pos(self, linenum):
        '''Reading initial positions of atoms from the line
        with index ``linenum``

        Args:
            linenum (int): the index of line 
                starting with 'position of ions in fractional'

        Returns:
            int, the index of ending line
        '''
        pos = []
        i = 0
        # positions of atoms
        while True:
            l = self._outlines[linenum+1+i].strip()
            if l == '':
                break
            pos.append(conv_string(l, float))
            i += 1
        self._initPos = pos
        return linenum + i

    def _read_init_latt(self, linenum):
        '''read the initial lattice vectors

        Args:
            linenum (int)

        Returns:
            int, the index of last line of the searched region
        '''
        latt = []
        for i in range(3):
            l = self._outlines[linenum+1+i]
            latt.append(conv_string(l, float, 0, 1, 2))
        self._initLatt = latt
        return linenum + 6

    def _read_weight(self, linenum):
        '''Read the weight of kpoints

        Args:
            linenum (int): the index of line starting with "Subroutine IBZKPT"

        Returns:
            int, the index of last line of the searched region
        '''
        nkpts = conv_string(self._outlines[linenum+3], int, 1)
        weight = []
        # for i in self._outlines[i+7:i+7+self.nirkp]:
        for i in range(nkpts):
            il = linenum + 7 + i
            l = self._outlines[il]
            weight.append(conv_string(l, int, -1))
        self._weight = weight
        return linenum + 3 + (nkpts+3) * 2 + 5

    @property
    def weight(self):
        return self._weight

    @property
    def nkpts(self):
        return self._nkpts

    def _read_kpts(self, linenum):
        '''Read the reciprocal coordinate of kpoints

        Args:
            linenum (int): the index of line 
                starting with "k-points in reciprocal lattice and weights"

        Returns:
            int, the index of last line of the searched region
        '''
        kpts = []
        weight = []
        i = 0
        while True:
            l = self._outlines[linenum+1+i].strip()
            if l == '':
                break
            kw = conv_string(l, float)
            kpts.append(kw[0:3])
            weight.append(kw[-1])
            i += 1
        self._kpoints = kpts
        if not hasattr(self, '_weight'):
            self._weight = weight
        return linenum + self._nkpts

    @property
    def kpoints(self):
        return self._kpoints

    def _read_dim_params(self, linenum):
        '''Read the dimension parameters of calculations

        Args:
            linenum (int): the index of line starting with 'Dimension of arrays'
        '''
        self._nkpts, self._nbands = conv_string(
            self._outlines[linenum+1], int, 3, -1)
        self._nedos = conv_string(self._outlines[linenum+2], int, 5)
        # FFT grids
        self._ngxyz = conv_string(self._outlines[linenum+6], int, 4, 7, -1)
        self._ngfxyz = conv_string(self._outlines[linenum+7], int, 3, 5, -1)
        # get the number of types and generate the symbols of atoms
        typeline = self._outlines[linenum+9].split('=')[-1]
        ntypes = conv_string(typeline, int)
        self._atoms = atoms_from_sym_nat(self._atomTypes, ntypes)
        return linenum + 15

    @property
    def atoms(self):
        return self._atoms

    @property
    def natoms(self):
        return len(self._atoms)

    @property
    def nbands(self):
        return self._nbands

    def _read_incar_params(self, linenum):
        '''Read import INCAR tags shown in OUTCAR
        '''
        tags = {}
        tagsOfInterest = {
            "ISPIN": (r"ISPIN", int, (2,)),
            "NELECT": (r"NELECT", float, (2,)),
            "ISMEAR": (r"ISMEAR", int, (2,)),
            "SIGMA": (r"ISMEAR", float, (5,)),
            "ENCUT": (r"ENCUT", float, (2,)),
        }
        i = linenum
        while True:
            l = self._outlines[i].strip()
            if l.startswith("SYSTEM"):
                tags["SYSTEM"] = trim_before(l, r"=").strip()
            if l.startswith("POSCAR"):
                self._poscarComm = trim_before(l, r"=").strip()
            for k, v in tagsOfInterest.items():
                if re.match(v[0], l):
                    tags[k] = conv_string(l, v[1], *v[2], strips=';')
            i += 1
            if l.startswith('-----'):
                break
        self._incar = Incar(**tags)
        return i

    @property
    def incar(self):
        return self._incar

    @property
    def nspins(self):
        return self._incar["ISPIN"]

    @property
    def encut(self):
        return self._incar["ENCUT"]

    def _read_npw(self, linenum):
        '''Read the number of plane waves at each k point
        '''
        self._npw = []
        for i in range(self.nkpts):
            l = self._outlines[linenum+i]
            self._npw.append(conv_string(l, int, -1))
        maxNpw, minNpw = conv_string(
            self._outlines[linenum+self.nkpts+1], int, -2, -1)
        # consistency check, related to the internal functionality
        assert maxNpw == max(self._npw)
        assert minNpw == min(self._npw)
        return linenum + self.nkpts + 3
    
    @property
    def maxNpw(self):
        return max(self._npw)

    def _init_iterations(self, linenum):
        '''initialize the regions of each ionic and electronic iterations

        Namely, it constructs a `_iterations` attribute as a 2d-array of int,
        each as the index of line in OUTCAR.
        The first dimension is the ionic step, and the second for SCF steps.
        The last two numbers in the second dimension marks the start and end
        of lines of the analysis of SCF loop.

        Note that the definition of iteration blocks by `_iterations` is only
        approximate, but adequate for further extraction of data.

        Args:
            linenum (int): the index of line starting with
                "First call to Ewald"

        Returns:
            int
        '''
        iters = []
        i = linenum + 1
        while True:
            l = self._outlines[i].strip().strip(' -')
            # electronic iteration
            if l.startswith('Iteration'):
                if conv_string(l, int, -1, strips=')') == 1:
                    scfiter = [i, ]
                else:
                    scfiter.append(i)
            # end of SCF of one ionic iteration
            if l.startswith('average (electrostatic)'):
                scfiter.append(i-1)
            # end of an ionic iteration
            if l.startswith('LOOP+'):
                scfiter.append(i+2)
                iters.append(scfiter)
            i += 1
            # break at timing, or EOF
            if l.startswith('General timing') or i == len(self._outlines):
                break
        self._iterations = iters
        self._nIonSteps = len(self._iterations)
        # ! need test for non-SCF cases
        self._nScfs = [len(x) - 2 for x in self._iterations]
        return i

    @property
    def nIonSteps(self):
        return self._nIonSteps

    def load_band(self, istep=None, kTrimBefore=None, kTrimAfter=None):
        '''Load the band structure at the end of a particular ion step.

        Note that wave projections are not available in OUTCAR. 
        For wave projections, use Vasprunxml class

        Args:
            istep (int): the index of ionic step.
                The last ionic step will be used if not specified.
            kTrimBefore (int): the kpoints before `kTrimBefore` will be trimed when parsing to BandStructure
            kTrimAfter (int): the kpoints after `kTrimAfter` will be trimed when parsing to BandStructure

        Returns:
            ``BandStructure`` object
        '''
        st, ed = self._get_ionic_data_region(istep=istep)
        eigen = []
        occ = []
        emult = {1: 2, 2: 1}[self.nspins]
        while st <= ed:
            l = self._outlines[st].strip()
            if l.startswith('E-fermi'):
                efermi = conv_string(l, float, 2)
                spanSpin = self.nkpts * \
                    (self.nbands + 3) + 2 * (self.nspins - 1)
                spanKpts = self.nbands + 3
                for spin in range(self.nspins):
                    # the index of line containing the first
                    # eigenvalue and occupation number of each spin
                    ln = st + 5 + 2 * (self.nspins - 1) + spanSpin * spin
                    eigenSpin = []
                    occSpin = []
                    for k in range(self.nkpts):
                        eoKpt = [conv_string(x, float, 1, 2)
                                 for x in self._outlines[ln+spanKpts*k:ln+spanKpts*k+self.nbands]]
                        eigenKpt, occKpt = tuple(np.transpose(eoKpt))
                    # the occupation numbers for ISPIN=1 have 2 as full occupancy
                    # divided by 2 for consistency with BandStructure class
                        for i, _v in enumerate(occKpt):
                            occKpt[i] /= emult
                        eigenSpin.append(eigenKpt)
                        occSpin.append(occKpt)
                    eigen.append(eigenSpin)
                    occ.append(occSpin)
                break
            st += 1

        stk = kTrimBefore
        if kTrimBefore is None:
            stk = 0
        edk = kTrimAfter
        if edk is None:
            edk = self.nkpts
        
        kvec = np.dot(self.kpoints, self.get_pos(istep=istep).b)
        return BandStructure(np.array(eigen, dtype=self._dtype)[:, stk:edk, :], 
                             np.array(occ, dtype=self._dtype)[:, stk:edk, :], 
                             weight=self._weight[stk:edk],
                             efermi=efermi, unit='ev',
                             kvec=kvec[stk:edk, :])

    def _get_ionic_data_region(self, istep=None):
        '''
        Return the region of data in self.outlines for a particular ionic step.

        Args:
            istep (int): The index of ionic step. 
                If not specified, the last step will be returned

        Returns:
            two int, the indices of starting and end line of data block of
            istep-th ionic interatioon.
        '''
        if istep is None:
            i = -1
        else:
            try:
                assert isinstance(istep, int)
                assert -self.nIonSteps < istep < self.nIonSteps
                i = istep
            except AssertionError:
                raise OutcarError(
                    "ionic step index {} overflows. (max {})".format(istep, self.nIonSteps-1))
        st = self._iterations[i][-2]
        ed = self._iterations[i][-1]
        return st, ed

    def get_pos(self, istep=None):
        '''
        Get the lattice structure for a particular ionic step.

        For initial structure, try ``initPoscar`` attribute

        Args:
            istep (int): the index of ionic step. 
                If not specified, the last ionic step will be used

        Returns:
            a ``Poscar`` object
        '''

        st, ed = self._get_ionic_data_region(istep)
        while st <= ed:
            l = self._outlines[st].strip()
            if l.startswith('direct lattice vectors'):
                latt = [conv_string(x, float, 0, 1, 2)
                        for x in self._outlines[st+1:st+4]]
            # the Cartisian coordinates
            if l.startswith('POSITION'):
                pos = [conv_string(x, float, 0, 1, 2)
                       for x in self._outlines[st+2:st+2+self.natoms]]
                break
            st += 1
        pc = Poscar(latt, self.atoms, pos, coordSys="C", unit="ang")
        # convert to direct coordinate
        pc.coordSys = "D"
        return pc

#     def get_total_force(self, iatom=None, istep=-1):
#         '''
#         Get the total-force information for a particular ionic step.

#         Parameters:
#             iatom: int
#                 The index of atom to output the total force
#             istep: int
#                 The index of ionic step. Default -1 to return the final step.
#                 istep=0 will give the same result as istep=1.

#         Returns:
#             total_force: numpy array
#                 The total force on all or one atom.
#                 If iatom is set as a non-negative value, shape(3). Otherwise shape(natoms, 3)
#         '''

#         if istep==0:
#             st, ed = self.__get_ionic_data_region(1, False)
#         else:
#             st, ed = self.__get_ionic_data_region(istep, False)

#         self.total_force = []

#         for i in range(st, ed):
#             #print(self.outlines[i])
#             if self._outlines[i].strip().startswith('POSITION'):
#                 for atom in range(self.natoms):
#                     self.total_force.append([float(x) for x in self._outlines[i+2+atom].split()[3:]])
#                 break

#         self.total_force = np.array(self.total_force)

#         if iatom is None:
#             return self.total_force
#         else:
#             try:
#                 assert int(iatom) < self.natoms
#             except ValueError:
#                 raise ValueError('invalid index of atom')
#             except AssertionError:
#                 raise ValueError('iatom out of range (%d)' % self.natoms)
#             else:
#                 return self.total_force[iatom]


def get_value(key, pathOut='OUTCAR'):
    '''
    Return the value in OUTCAR by keyword.

    Args:
        key (str): the keyword of value to extract. Case insensitive.
        Accepted are
            "mnpw": maximum number of plane waves
            "nbands": NBANDS
            "enuct": ENCUT
            "nkpts": number of kpoints in calculation
            "efermi": Fermi-energy
        pathOut (str): the path of OUTCAR file
    '''
    k = key.lower()
    oc = Outcar(pathOut=pathOut)
    kv = {
        "mnpw": oc.maxNpw,
        "nbands": oc.nbands,
        "encut": oc.encut,
        "nkpts": oc.nkpts,
        "efermi": oc.efermi
    }
    return kv.get(k, None)
