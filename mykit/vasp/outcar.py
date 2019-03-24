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
from mykit.core.utils import conv_string
from mykit.vasp.incar import Incar
from mykit.vasp.poscar import Poscar


class OutcarError(Exception):
    pass


class Outcar(Verbose):

    def __init__(self, pathOut='OUTCAR'):

        # read in OUTCAR
        with open(pathOut, 'r') as f:
            self._outlines = f.readlines()

        # self.eigenene = None
        # self.occ = None
        self._check_finished()
        self._read_all()
        # Build initial POSCAR
        self._initPoscar = Poscar(self._initLatt, self.atoms, self._initPos,
                                  unit="ang", coordSys="D")
        # self._divide_ion_iteration()

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
            r'k-points in reciprocal': self._read_kpts,
            r'position of ions in fractional': self._read_init_pos,
            r'Dimension of arrays': self._read_dim_params,
            r'SYSTEM =': self._read_incar_params,
            r'k-point  1 :': self._read_npw,
            r'total amount of memory': self._init_iterations,
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
        self._nkpts = len(weight)
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
                starting with "k-points in reciprocal lattice"

        Returns:
            int, the index of last line of the searched region
        '''
        kpts = []
        for i in range(self.nkpts):
            l = self._outlines[linenum+1+i]
            kpts.append(conv_string(l, float, 0, 1, 2))
        self._kpoints = kpts
        return linenum + self.nkpts

    @property
    def kpoints(self):
        return self._kpoints

    def _read_dim_params(self, linenum):
        '''Read the dimension parameters of calculations

        Args:
            linenum (int): the index of line starting with 'Dimension of arrays'
        '''
        nkpts, self._nbands = conv_string(
            self._outlines[linenum+1], int, 3, -1)
        assert nkpts == self.nkpts
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

    def _read_incar_params(self, linenum):
        '''Read import INCAR tags shown in OUTCAR
        '''
        return linenum

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

    def _init_iterations(self, linenum):
        '''initialize the regions of each ionic and electronic iterations
        '''
        return linenum

#         # the parameters are saved in the part before the iterations start
#         for i in range(len(self._outlines[:self.iterations[0][0]])):
#         #for i in range(len(self.outlines)):
#             line = self._outlines[i].strip()


#         # NBANDS, NGX/Y/Z
#             if line.startswith('Dimension of arrays'):
#                 self.nbands = int(self._outlines[i+1].split()[-1])
#                 self.ngxyz  = [int(x) for x in self._outlines[i+6].split()[4:10:3]]
#                 self.ngfxyz  = [int(x) for x in self._outlines[i+7].split()[3:7:2]]
#                 self.natoms = int(self._outlines[i+2].split()[-1])
#         # Maximum number of planewaves
#             if line.startswith('maximum number of plane-waves'):
#                 self.mplw = conv_string(line, int, -1)
#         # ISPIN
#             if line.startswith('ISPIN'):
#                 self.ispin = conv_string(line, int, 2)
#         # NELECT: number of electrons
#             if line.startswith('NELECT'):
#                 self.nelect = conv_string(line, int, 2)
#         # ISMEAR and SIGMA: broadening information
#             if line.startswith('ISMEAR'):
#                 self.ismear = conv_string(line.replace(';',''), int, 2)
#                 self.sigma  = conv_string(line.replace(';',''), float, 5)

#         # lattice constants and inner coordinates
#         self.lattice, self.innerpos = self.get_pos(ionstep=0)

#     def _divide_ion_iteration(self):
#         '''
#         Check the line number of each ionic step
#         '''
#         ln_iteration = []
#         i_ionstep = []
#         self.flag_static = True

#         for i in range(len(self._outlines)):
#             words = self._outlines[i].split()
#             if len(words) == 5:
#                 if words[1] == 'Iteration':
#                     ln_iteration.append(i)
#                     i_ionstep.append(int(words[2][:-1]))

#         self.nisteps = i_ionstep[-1]
#         if self.nisteps == 1:
#             self.iterations = ln_iteration
#         else:
#             self.flag_static = False
#             self.iterations = [[] for i in range(self.nisteps)]
#             for i in range(len(ln_iteration)):
#                 self.iterations[i_ionstep[i]-1].append(ln_iteration[i])


#     def __get_ionic_data_region(self, ionstep, all_elesteps=False):
#         '''
#         Return the region of data in self.outlines for a particular ionic step.

#         Parameters:
#             ionstep: int
#                 The index of ionic step. Default -1 to return the information regarding the final results.
#                 0 for the initial data.
#             all_elesteps: bool
#                 flag for controlling whether the whole range of electronic step or only the last step is included.
#         '''

#         if all_elesteps:
#             end_ele = 0
#         else:
#             end_ele = -1

#         if ionstep == 0:
#             st_line = 0
#             ed_line = self.iterations[0][0]
#         elif ionstep == -1 or ionstep == self.nisteps:
#             st_line = self.iterations[-1][end_ele]
#             ed_line = len(self._outlines)-1
#         elif ionstep < -1 or ionstep > self.nisteps:
#             raise ValueError("The requested ionic step is too large. Maximum %d" % self.nisteps)
#         else:
#             st_line = self.iterations[ionstep-1][end_ele]
#             ed_line = self.iterations[ionstep][0]

#         return st_line, ed_line


#     def get_pos(self, ionstep=-1):
#         '''
#         Get the initial lattice structure for a particular ionic step. Default return the last structure.

#         Parameters:
#             ionstep: int
#                 The index of ionic step. Default -1 to return the final structure.
#                 0 for the initial structure, i.e. POSCAR
#         '''

#         st_line, ed_line = self.__get_ionic_data_region(ionstep, False)
#         innerpos = []

#         for i in range(st_line, ed_line):
#             #print(self.outlines[i])
#             if self._outlines[i].strip().startswith('direct'):
#                 lattice = np.array([
#                                 [float(x) for x in self._outlines[i+1].split()[:3]],
#                                 [float(x) for x in self._outlines[i+2].split()[:3]],
#                                 [float(x) for x in self._outlines[i+3].split()[:3]],
#                                    ])
#             if ionstep == 0:
#                 if self._outlines[i].strip().startswith('position of ions in'):
#                     for atom in range(self.natoms):
#                         innerpos.append([float(x) for x in self._outlines[i+1+atom].split()])
#                     break
#             else:
#                 if self._outlines[i].strip().startswith('POSITION'):
#                     for atom in range(self.natoms):
#                         innerpos.append([float(x) for x in self._outlines[i+2+atom].split()[:3]])
#                     break

#         innerpos = np.array(innerpos)
#         # convert cartisian to direct coordinates
#         if ionstep != 0:
#             innerpos = np.dot(innerpos, np.linalg.inv(lattice))
#             # the inner coordinates are rounded, which may cause numerical error
#             innerpos = np.round(innerpos,5)
#         return lattice, innerpos


#     def get_total_force(self, iatom=None, ionstep=-1):
#         '''
#         Get the total-force information for a particular ionic step.

#         Parameters:
#             iatom: int
#                 The index of atom to output the total force
#             ionstep: int
#                 The index of ionic step. Default -1 to return the final step.
#                 ionstep=0 will give the same result as ionstep=1.

#         Returns:
#             total_force: numpy array
#                 The total force on all or one atom.
#                 If iatom is set as a non-negative value, shape(3). Otherwise shape(natoms, 3)
#         '''

#         if ionstep==0:
#             st_line, ed_line = self.__get_ionic_data_region(1, False)
#         else:
#             st_line, ed_line = self.__get_ionic_data_region(ionstep, False)

#         self.total_force = []

#         for i in range(st_line, ed_line):
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


#     def get_band_structure(self, ionstep=-1):
#         '''
#         Get the band structure for a particular ionic step.

#         Parameters:
#             ionstep: int
#                 The index of ionic step. -1 to return the band structure of final structure
#                 0 for that of the initial structure
#         '''
#         st_line, ed_line = self.__get_ionic_data_region(ionstep, False)
#         # find the starting line for the band structure
#         for i in range(st_line,ed_line):
#             line = self._outlines[i].strip()
#             if line.startswith('E-fermi'):
#                 ln_e_fermi = i
#                 break

#         nbands = self.nbands
#         ispin = self.ispin
#         nirkp = self.nirkp

#         kpt_block = nbands + 3
#         spin_block = kpt_block * nirkp + 2

#         self.efermi = conv_string(self._outlines[ln_e_fermi], float, 2)
#         self.eigenene = []
#         self.occ = []
#         vb = self.nelect / 2 - 1
#         cb = self.nelect / 2
#         self.vbm = -10000.0
#         self.cbm =  10000.0
#         self.vbm_loc =  [0, 0]
#         self.cbm_loc =  [0, 0]
#         self.Eg_dir = 10000.0
#         self.dir_loc =  [0, 0]

#         if ispin == 1:
#             ln_spin_1_k_1 = ln_e_fermi + 3
#         else:
#             ln_spin_1_k_1 = ln_e_fermi + 5 # extra "spin component" and empty lines

#         for spin in range(ispin):
#             e_spin = []
#             occ_spin = []
#             for kpt in range(nirkp):
#                 e_kpt = []
#                 occ_kpt = []
#                 band_st_line = ln_spin_1_k_1 + spin*spin_block + kpt*kpt_block + 2
#                 for line in self._outlines[band_st_line:band_st_line+nbands]:
#                     words = line.split()
#                     e_kpt.append(float(words[1]))
#                     occ_kpt.append(float(words[2]))

#                 Eg_dir = e_kpt[cb] - e_kpt[vb]
#                 if Eg_dir < self.Eg_dir:
#                     self.Eg_dir  = Eg_dir
#                     self.dir_loc = [spin, kpt]

#                 if e_kpt[vb] > self.vbm:
#                     self.vbm     = e_kpt[vb]
#                     self.vbm_loc = [spin, kpt]
#                 if e_kpt[cb] < self.cbm:
#                     self.cbm     = e_kpt[cb]
#                     self.cbm_loc = [spin, kpt]

#                 e_spin.append(e_kpt)
#                 occ_spin.append(e_kpt)
#             self.eigenene.append(e_spin)
#             self.occ.append(occ_spin)

#         # the fundamental band gap
#         self.Eg = self.cbm - self.vbm
#         self.Eg_vbm = self.eigenene[self.vbm_loc[0]][self.vbm_loc[1]][cb] - self.eigenene[self.vbm_loc[0]][self.vbm_loc[1]][vb]
#         self.Eg_cbm = self.eigenene[self.cbm_loc[0]][self.cbm_loc[1]][cb] - self.eigenene[self.cbm_loc[0]][self.cbm_loc[1]][vb]
#         self.eigenene = np.array(self.eigenene)
#         self.occ = np.array(self.occ)

#         return self.eigenene, self.occ


# def read_fermi(outcar='OUTCAR'):

#     fermi_level = 0.0

#     if os.path.exists(outcar):
#         try:
#             fermi_level = float(sp.check_output("awk '/E-fermi/ {print $3}' %s" % outcar, shell=True))
#             print("Fermi level found: %6.4f" % fermi_level)
#         except ValueError:
#             print("WARNING: error in reading E-fermi from %s. Fermi level is set to 0." % outcar)
#     else:
#         print("WARNING: %s is not found. Fermi level is set to 0." % outcar)

#     return fermi_level


# def vasp_anal_get_outcar(keyin,index=-1,outcar='OUTCAR'):
#     '''
#     return the value of key in outcar.
#     '''
#     key = keyin.lower()
# # maximum number of plane wave, i.e. the maximum size of representation matrix
#     if key=='mnpw':
#         mnpw = sp.check_output("awk '/maximum number of/ {print $5}' %s | tail -1" % outcar,shell=True)
#         return int(mnpw)
# # NBANDS
#     if key in ['nb', 'nbands']:
#         nb = sp.check_output("awk '/NBANDS/ {print $15}' %s | head -1" % outcar,shell=True)
#         return nb
# # ENCUT
#     if key in ['encut']:
#         encut = sp.check_output("awk '/ENCUT/ {print $3}' %s | head -1" % outcar,shell=True)
#         encut = int(float(encut))
#         return encut
# # converged G.S. energy
#     if key in ['ene', 'energy']:
#         ene = sp.check_output("awk '/without/ {print $7}' %s | tail -1" % outcar,shell=True)
#         ene = float(ene)
#         return ene
# # total number of irreducible k-points
#     if key in ['nkp']:
#         nkp = sp.check_output("awk '/NKPTS/ {print $4}' %s | head -1" % outcar,shell=True)
#         nkp = int(nkp)
#         return nkp
#     if key in ['efermi', 'e-fermi', 'fermi']:
#         efermi = float(sp.check_output("awk '/E-fermi/ {print $3}' %s | tail -1" % outcar, shell=True))
#         return efermi
