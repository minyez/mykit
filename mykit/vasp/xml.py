# coding = utf-8
'''
'''

import re
from copy import deepcopy

import numpy as np

from mykit.core.bandstructure import BandStructure
from mykit.core.cell import sym_nat_from_atoms
from mykit.core.log import Verbose
from mykit.core.numeric import Prec
from mykit.core.utils import conv_string, get_first_last_line, get_str_indices
from mykit.vasp.incar import Incar
from mykit.vasp.poscar import Poscar

try:
    from lxml import etree
except ModuleNotFoundError:
    from xml.etree import ElementTree as etree


class VasprunxmlError(Exception):
    pass


class Vasprunxml(Verbose, Prec):
    '''Class to read and analyse the data from the vasprun.xml

    The shape of import properties:
        eigen (array): shape (nspins, nibzkpt, nbands)
        occ (array): shape (nspins, nibzkpt, nbands)
        dos (array): shape (nspins, nedos)
        totalDos: alias of dos

    Optional properties, depending on the input:
        pDos (array): shape (natoms, nspins, nedos, nprojs)
        pWave (array): shape (nspins, nibzkpt, nbands, natoms, nprojs)
    '''

    def __init__(self, pathXml='vasprun.xml', *args):
        flag = _check_vasprunxml_status(pathXml)
        if flag is None:
            raise VasprunxmlError("Parsed file seems not a vasprun.xml")
        elif not flag:
            self.print_warn(
                "the calculation is not finished. Try to read it...")
        tree = etree.parse(pathXml)
        self.__root = tree.getroot()
        self._init_sections()
        self._read_atoms()
        self._initPoscar = self._read_pos(self._secInitStruct,
                                          comment="vasprun.xml initialpos")
        self._finalPoscar = self._read_pos(self._secFinalStruct,
                                           comment="vasprun.xml finalpos")
        self._read_incar()
        self._read_params()
        self._read_klist()
        self._read_eigen_occ()
        self._hasProjected = False
        self._read_dos()
        if self._hasProjected:
            self._read_projected()

    def _init_sections(self):
        '''Initialize the XML roots of each sections
        '''
        # Manual parameters by INCAR
        self._secIncar = self.__root.find('incar')
        # Initial and final structures
        self._secInitStruct = self.__root.find(
            './/structure[@name="initialpos"]')
        self._secFinalStruct = self.__root.find(
            './/structure[@name="finalpos"]')
        # auto-generated parameters
        self._secPara = self.__root.find('parameters')
        # all ionic iterations
        self._secCalcs = self.__root.findall('calculation')
        # the last ionic calculation has the final eigenvalues
        self._secCalcLast = self._secCalcs[-1]
        self._secDos = self._secCalcLast.find('dos')
        self._secAtoms = self.__root.find('atominfo')
        self._secKpoints = self.__root.find('kpoints')

    def _read_pos(self, secStructure, comment=None):
        '''Return Poscar from a vasprun.xml structure tag section

        In vasprun.xml, the positions are always stored in direct system.

        Returns:
            Poscar instance
        '''
        lattVArray = secStructure.find(
            'crystal').find('.//varray[@name="basis"]')
        latt = [list(map(float, x.text.split())) for x in lattVArray]
        posVArray = secStructure.find('varray')
        pos = [list(map(float, x.text.split())) for x in posVArray]
        return Poscar(latt, self.atoms, pos, coordSys="D", comment=comment)

    def _read_incar(self):
        '''Read manual parameters from incar section
        '''
        tags = {}
        for i in self._secIncar:
            incarLine = i.attrib["name"] + "=" + i.text
            tags.update(Incar.analyze_incar_line(incarLine))
        self.incar = Incar(**tags)

    @property
    def nelect(self):
        return self._nelect

    @property
    def nspins(self):
        return self._nspins

    @property
    def nbands(self):
        return self._nbands

    def _read_params(self):
        '''Read auto-generated parameters
        '''
        # ISPIN: root->parameter->separator name='electronic'->separator name='elecronic spin'->[0]
        # NBANDS: root->parameter->separator name='electronic'->i name='NBANDS'
        paramE = self._secPara.find('.//separator[@name="electronic"]')
        spin = paramE.find('.//separator[@name="electronic spin"]')
        self._nelect = float(paramE.find('.//i[@name="NELECT"]').text)
        self._nspins = int(spin[0].text)
        self._nbands = int(paramE.find('.//i[@name="NBANDS"]').text)

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

    def get_atom_index(self, atomType=None):
        '''Get the indices of atoms of ``atomType``.

        If atomType is not specified, all indices will be returned.
        If atomType is int, the indices of atoms with atomType-th type will be returned

        Args:
            atomType (str or int): the identifier of atomic type
        '''
        at = atomType
        if at is None:
            return list(range(self.natoms))
        if isinstance(at, int):
            pass
        elif isinstance(at, str):
            try:
                at = self.atomTypes.index(at)
            except ValueError:
                raise VasprunxmlError("Atom type not found: {}".format(at))
        else:
            raise VasprunxmlError("Atomtype should be either int or string.")
        at = self.atomTypes[at]
        return get_str_indices(self._atoms, at)

    @property
    def eigen(self):
        return self._eigen

    @property
    def occ(self):
        return self._occ

    def _read_eigen_occ(self):
        # self._eigen = np.zeros([self.nspins,self._nibzkpt,self.nbands])
        self._eigen = []
        self._occ = []
        eigenSet = self._secCalcLast.find(
            'eigenvalues').find('array').find('set')
        for spin in range(self.nspins):
            eigenSp = []
            occSp = []
            for kp in range(self._nibzkpt):
                # eo = [list(map(float, x.text.split())) for x in eigenSet[spin][kp]]
                eo = [conv_string(x.text, float) for x in eigenSet[spin][kp]]
                eo = np.array(eo)
                eigenSp.append(eo[:, 0])
                occSp.append(eo[:, 1])
            self._eigen.append(eigenSp)
            self._occ.append(occSp)
        # self._eigen = np.array(self._eigen)
        # self._occ = np.array(self._occ)

    @property
    def pDos(self):
        if hasattr(self, '_pDos'):
            return self._pDos
        return None

    @property
    def dosGrid(self):
        if hasattr(self, '_dosGrid'):
            return self._dosGrid
        return None

    @property
    def totalDos(self):
        if hasattr(self, '_totalDos'):
            return self._totalDos
        return None

    @property
    def dos(self):
        return self.totalDos

    @property
    def nprojs(self):
        if hasattr(self, '_projs'):
            return len(self._projs)
        return 0

    @property
    def projs(self):
        if hasattr(self, '_projs'):
            return self._projs
        return None

    @property
    def efermi(self):
        return self._efermi

    def _read_dos(self):
        '''
        TODO:
            read in the DOS
        '''
        sd = self._secDos
        self._NEDOS = 0
        if sd is None:
            self._efermi = None
        else:
            self._efermi = float(sd[0].text)
            totDosSet = sd.find('total').find('array').find('set')
            self._dosGrid = None
            self._totalDos = []
            self._totalDosInteg = []
            for spin in range(self.nspins):
                v = [conv_string(x.text, float) for x in totDosSet[spin]]
                v = np.transpose(v)
                if spin == 0:
                    self._dosGrid = v[0]
                    self._NEDOS = len(v[0])
                self._totalDos.append(v[1])
                self._totalDosInteg.append(v[2])
            # attemp to read projected DOS
            pDosRoot = sd.find('partial')
            if pDosRoot != None:
                self._hasProjected = True
            if self._hasProjected:
                pDosArray = pDosRoot.find('array')
                # projector strings: root->dos->partial->array->field[1:]
                # the first is energy label
                self._projs = [i.text.strip()
                               for i in pDosArray.findall('field')[1:]]
                pDosSet = pDosArray.find('set')
                _pDos = []
                for ia in range(self.natoms):
                    _pDos.append([])
                    for spin in range(self.nspins):
                        v = [conv_string(x.text, float)[1:]
                             for x in pDosSet[ia][spin]]
                        _pDos[-1].append(v)
                # original shape (natoms, nspins, nedos, nprojs)
                # swap axis for  (nspins, nedos, natoms, nprojs)
                self._pDos = np.swapaxes(np.swapaxes(_pDos, 0, 1), 1, 2)

    @property
    def pWave(self):
        if hasattr(self, '_pWave'):
            return self._pWave
        return None

    def _read_projected(self):
        projSec = self._secCalcLast.find('projected')
        # pWave data: root->calculation->projected->array->set
        pWaveSet = projSec.find('array').find(
            'set')  # contains ISPIN dataset(s)

        # initialize the pWave data
        # self._pWave = np.zeros([self.nspins,self._nibzkpt,self.nbands,self.natoms,len(self._projs)])
        _pWave = []
        # shape: (nspin, nibzkpt, nbands, natoms, nprojs)
        for spin in range(self.nspins):
            _pWave.append([])
            _listSp = _pWave[-1]
            for kp in range(self._nibzkpt):
                _listSp.append([])
                _listKp = _listSp[-1]
                for band in range(self.nbands):
                    _listKp.append([conv_string(at.text, float)
                                    for at in pWaveSet[spin][kp][band]])
        self._pWave = _pWave

    def load_band(self, kTrimBefore=None, kTrimAfter=None):
        '''Return a BandStructure instance.

        Note:
            `BandStructure.nelect` attribute may differ from `Vasprunxml.nelect`, 
            if the smearing is large when compared with the band gap.

        Args:
            kTrimBefore (int): the kpoints before `kTrimBefore` will be trimed when parsing to BandStructure
            kTrimAfter (int): the kpoints after `kTrimAfter` will be trimed when parsing to BandStructure
        '''
        stk = kTrimBefore
        if kTrimBefore is None:
            stk = 0
        edk = kTrimAfter
        if edk is None:
            edk = self.nibzkpt
        projected = None
        if self.pWave != None:
            projected = {
                "atoms": self._atoms,
                "projs": self.projs,
                "pWave": np.array(self.pWave, dtype=self._dtype)[:, stk:edk, :, :, :]}
        bs = BandStructure(
            np.array(self._eigen, dtype=self._dtype)[:, stk:edk, :],
            np.array(self._occ, dtype=self._dtype)[:, stk:edk, :],
            np.array(self._weight, dtype=self._dtype)[stk:edk],
            efermi=self._efermi, projected=projected,
            kvec=np.array(self.kvec, dtype=self._dtype)[stk:edk, :],
        )
        return bs

    @property
    def nibzkpt(self):
        return self._nibzkpt

    @property
    def nkpt(self):
        return self._nkpt

    @property
    def kvec(self):
        '''list. the vectors of all kpoints

        Reciprocal basis of final lattice is used.
        '''
        return np.dot(self._ibzkpt, np.transpose(self._finalPoscar.b))

    @property
    def kpoints(self):
        return self._ibzkpt

    @property
    def weight(self):
        return self._weight

    @property
    def kptsWeight(self):
        return [kpt + [w, ] for kpt, w in zip(self._ibzkpt, self._weight)]

    @property
    def kdense(self):
        if hasattr(self, "_kdense"):
            return self._kdense
        return None

    @property
    def kmode(self):
        if hasattr(self, "_kmode"):
            return self._kmode
        return None

    @property
    def kdiv(self):
        if hasattr(self, "_kdiv"):
            return self._kdiv
        return None

    def _read_klist(self):
        #       klist_index is 1 if auto generator is used
        #       or 0 if mannually included
        self._nkpt = 0
        gen = self._secKpoints.find('generation')
        if gen is None:
            ki = 0
        else:
            ki = 1
            kmode = gen.attrib["param"]
            # kdiv is only applicable for G and M
            if kmode in ["Gamma", "Monkhorst-Pack"]:
                self._kdiv = conv_string(gen[0].text, int)
                self._nkpt = np.prod(self._kdiv)
            elif kmode == "listgenerated":
                self._kdense = int(gen[0].text)
            elif kmode == "Auto":
                self._kdense = int(gen[0].text)
                self._kdiv = conv_string(gen[1].text, int)
                self._nkpt = np.prod(self._kdiv)
            else:
                raise VasprunxmlError("Unknown kmode")
            self._kmode = kmode[0].upper()
        self._nibzkpt = len(self._secKpoints[ki])
        if self._nkpt == 0:
            self._nkpt = deepcopy(self._nibzkpt)
        self._ibzkpt = [conv_string(kp.text, float)
                        for kp in self._secKpoints[ki]]
        self._weight = [int(np.rint(float(kp.text) * self._nkpt))
                        for kp in self._secKpoints[ki+1]]

    # def pwave_index(self,lcomponent):
    #     if self.projs == ():
    #         raise ValueError("XML does not contain partial wave information.")
    #     index_l = []
    #     # total wave
    #     if lcomponent == 't':
    #         index_l = list(range(len(self._projs)))
    #     for x in self._projs:
    #         if x.startswith(lcomponent):
    #             index_l.append(self._projs.index(x))
    #     return index_l

#     def sum_atom_l_comp(self, spin, band, kp, at_index, pw_index):
#         '''Sum the component of the atoms in at_index and partial waves in pw_index
#         '''
#         weigh = 0
#         for at in at_index:
#             for pw in pw_index:
# #                weigh += self._pWave[spin][band][kp][at][pw]
#                 weigh += self._pWave[spin][band][kp][at][pw]
#         return weigh

#     def get_gap(self):
#         # return the gap from eigenvalue information
#         ecbm = 100000.0
#         evbm = -100000.0
#         vbm = self.nelect/2
#         cbm = self.nelect/2 + 1
#         for spin in range(self.nspins):
#             for kp in range(self._nibzkpt):
# #                print self.eigen[spin,kp,vbm-1],self.eigen[spin,kp,cbm-1]
#                 if (self.eigen[spin,kp,vbm-1]> evbm):
#                     evbm = self.eigen[spin,kp,vbm-1]
#                 if (self.eigen[spin,kp,cbm-1]< ecbm):
#                     ecbm = self.eigen[spin,kp,cbm-1]
#         gap = ecbm - evbm
#         return gap


def _check_vasprunxml_status(pathXml):
    '''Check if the calculation of vasprun.xml ``pathXml`` has finished itself.

    Args:
        pathXml (str): the path of vasprun.xml to check

    Returns:
        None, if pathXml does not exist or is not an vasp XML file
        True if the calculation finished, otherwise False
    '''
    try:
        first, last = get_first_last_line(pathXml)
    except FileNotFoundError:
        return None
    if not re.match(r"<\?xml version=\"[\d.]+\" encoding=\"[\w-]+\"\?>",
                    first.strip()):
        return None
    if last.strip() == "</modeling>":
        return True
    else:
        return False
