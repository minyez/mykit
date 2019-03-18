# coding = utf-8
'''
'''

from copy import deepcopy

import numpy as np

from mykit.core.cell import get_sym_index, sym_nat_from_atoms
from mykit.core.log import Verbose
from mykit.core.numeric import Prec
from mykit.vasp.incar import Incar
from mykit.vasp.poscar import Poscar

try:
    from lxml import etree
except ModuleNotFoundError:
    import elementtree.ElementTree as etree


class VasprunxmlError(Exception):
    pass


class Vasprunxml(Verbose, Prec):
    '''Class to read and analyse the data from the vasprun.xml
    '''

    def __init__(self, xmlFile='vasprun.xml', *args):
        tree = etree.parse(xmlFile)
        self.__root = tree.getroot()
        self._init_sections()
        self._read_atoms()
        self._read_params()
        self._read_klist()
        self._read_eigen_occ()
        try:
            if 'pw' in args:
                self._read_pwdata()
        except AttributeError:
            self.print_warn("no PDOS data found")

    def _init_sections(self):
        self._secInitStruct = self.__root.find('.//structure[@name="initpos"]')
        self._secFinalStruct = self.__root.find('.//structure[@name="finalpos"]')
        self._secIncar = self.__root.find('incar')
        # parameters
        self._secPara = self.__root.find('parameters')
        # all ionic calculations
        self._secCalcs = self.__root.findall('calculation')
        # the last ionic calculation which has the final eigenvalues
        self._secCalcLast = self._secCalcs[-1]
        self._secAtoms = self.__root.find('atominfo')
        self._secKpoints = self.__root.find('kpoints')

    def _read_params(self):
        # ISPIN: root->parameter->separator name='electronic'->separator name='elecronic spin'->[0]
        # NBANDS: root->parameter->separator name='electronic'->i name='NBANDS'
        para_e = self._secPara.find('.//separator[@name="electronic"]')
        spin = para_e.find('.//separator[@name="electronic spin"]')
        self.nelec = int(float(para_e.find('.//i[@name="NELECT"]').text))
        self.ispin = int(spin[0].text)
        self.nbands = int(para_e.find('.//i[@name="NBANDS"]').text)

    def _read_pwdata(self):
        # partial wave strings: root->dos->partial->array->field[1:] # the first is energy label
        partial_array = self._secCalcLast.find('dos').find('partial').find('array')
        self.str_pwaves = [pwave.text.strip() for pwave in partial_array.findall('field')[1:]]
        # pdos data: root->calculation->projected->array
        proj =  self._secCalcLast.find('projected')
        pwav_dataset = proj.findall('array')[0][-1] # contains ISPIN dataset(s)

        # initialize the pwav data
        self.pwdata = np.zeros([self.ispin,self._nibzkpt,self.nbands,self.natoms,len(self.str_pwaves)])

        for spin in range(self.ispin):
            for kp in range(self._nibzkpt):
                for band in range(self.nbands):
                    for atom in range(self.natoms):
                        # data = np.array([float(x) for x in pwav_dataset[spin][kp][band][atom].text.split()])
                        self.pwdata[spin, kp, band, atom, :] = np.array([float(x) for x in pwav_dataset[spin][kp][band][atom].text.split()])

    @property
    def natoms(self):
        return len(self._atoms)
    @property
    def atoms(self):
        return self._atoms
    @property
    def typeMapping(self):
        '''Map index (int) to atom type (str)
        '''
        _ats = self.atomTypes
        _dict = {}
        for i, _at in enumerate(_ats):
            _dict.update({i: _at})
        return _dict
    @property
    def atomTypes(self):
        return self._atomTypes
    @property
    def natomsPerType(self):
        return self._natomsPerType
    @property
    def ntypes(self):
        return len(self._atomTypes)

    def _read_atoms(self):
        self._atoms = []
        atomsSet = self._secAtoms.find('.//array[@name="atoms"]').find('set')
        for rc in atomsSet:
            self._atoms.append(rc[0].text.strip())
        self._atomTypes, self._natomsPerType = sym_nat_from_atoms(self._atoms)

    def atoms_index(self, atomtype=None):
        if atomtype is None:
            return list(range(self.natoms))
        _at = None
        if isinstance(atomtype, int):
            try:
                assert 0 <= atomtype < self.ntypes
            except AssertionError:
                raise VasprunxmlError("Invalid atomtype index: {}".format(atomtype))
            _at = atomtype
        elif isinstance(atomtype, str):
            try:
                _at = self.atomTypes.index(atomtype)
            except ValueError:
                raise VasprunxmlError("Atom type not found: {}".format(atomtype))
        else:
            raise VasprunxmlError("Atomtype should be either int or string.")
        return get_sym_index(self._atoms, _at)

    @property
    def eigen(self):
        return self._eigen
    @property
    def occ(self):
        return self._occ

    def _read_eigen_occ(self):
        # self._eigen = np.zeros([self.ispin,self._nibzkpt,self.nbands])
        self._eigen = []
        self._occ = []
        eigenSet= self._secCalcLast.find('eigenvalues').find('array').find('set')
        for spin in range(self.ispin):
            eigenSp = []
            occSp = []
            for kp in range(self._nibzkpt):
                eo = [list(map(float, x.text.split())) for x in eigenSet[spin][kp]]
                eo = np.array(eo)
                eigenSp.append(eo[:, 0])
                occSp.append(eo[:, 1])
            self._eigen.append(eigenSp)
            self._occ.append(occSp)
        self._eigen = np.array(self._eigen)
        self._occ = np.array(self._occ)

    def _read_klist(self):
#       klist_index is 1 if auto generator is used
#       or 0 if mannually included
        self._nkpt = 0
        gen = self._secKpoints.find('generation')
        if gen is None:
            ki = 0
        else:
            ki = 1
            self._nkpt = np.prod(list(map(int, gen[0].text.split())))
        self._nibzkpt = len(self._secKpoints[ki])
        if self._nkpt == 0:
            self._nkpt = deepcopy(self._nibzkpt)
        # self._ibzkpt = [[float(kvec) for kvec in kp.text.split()] for kp in self._secKpoints[ki]]
        self._ibzkpt = [list(map(float, kp.text.split())) for kp in self._secKpoints[ki]]
        self._weight = [int(np.rint(float(kp.text) * self._nkpt)) for kp in self._secKpoints[ki+1]]

    def pwave_index(self,lcomponent):
        index_l = []
        # total wave
        if lcomponent == 't':
            index_l = list(range(len(self.str_pwaves)))
        for x in self.str_pwaves:
            if x.startswith(lcomponent):
                index_l.append(self.str_pwaves.index(x))
        return index_l

# sum the component of the atoms in at_index and partial waves in pw_index
# can be more pythonic
    def sum_atom_l_comp(self, spin, band, kp, at_index, pw_index):
        weigh = 0
        for at in at_index:
            for pw in pw_index:
#                weigh += self.pwdata[spin][band][kp][at][pw]
                weigh += self.pwdata[spin][band][kp][at][pw]
        return weigh

#     def get_gap(self):
#         # return the gap from eigenvalue information
#         ecbm = 100000.0
#         evbm = -100000.0
#         vbm = self.nelec/2
#         cbm = self.nelec/2 + 1
#         for spin in range(self.ispin):
#             for kp in range(self._nibzkpt):
# #                print self.eigen[spin,kp,vbm-1],self.eigen[spin,kp,cbm-1]
#                 if (self.eigen[spin,kp,vbm-1]> evbm):
#                     evbm = self.eigen[spin,kp,vbm-1]
#                 if (self.eigen[spin,kp,cbm-1]< ecbm):
#                     ecbm = self.eigen[spin,kp,cbm-1]
#         gap = ecbm - evbm
#         return gap
