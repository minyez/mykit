# coding = utf-8
'''Module that defines classes for VASP output files
'''
import os
import subprocess as sp

import numpy as np

# from mykit.core.lattice import lattice
from mykit.core.log import Verbose
from mykit.core.utils import conv_string


class OutcarError(Exception):
    pass


class Outcar(Verbose):
    
    def __init__(self, OutFile='OUTCAR', verbose=True):
        
        self.filename = OutFile
        self.eigenene = None
        self.occ = None
        self.verbose = verbose

        self.__print(" Reading OUTCAR: %s" % OutFile)
        with open(OutFile,'r') as f:
            self.outlines = f.readlines()
        self.__divide_ion_iteration()
        self.__init_calc_params()
        self.__check_finished()
        if not self.flag_finish:
            self.print_warn("Task has not finished.", level=0)

    def __print(self, pstr):
        self.print_log(pstr, self.verbose)

    def __init_calc_params(self):
        '''
        Initialize the parameters of calculations
        '''
        # the parameters are saved in the part before the iterations start
        for i in range(len(self.outlines[:self.iterations[0][0]])):
        #for i in range(len(self.outlines)):
            line = self.outlines[i].strip()

        # irreducible k-points
            if line.startswith('Subroutine IBZKPT returns following result'):
                self.nirkp = int(self.outlines[i+3].split()[1])
                self.irkp = []
                self.wirkp = []
                sum_weight = 0.0
                for k_str in self.outlines[i+7:i+7+self.nirkp]:
                    k_str_list = k_str.split()
                    self.irkp.append([float(x) for x in k_str_list[:3]])
                    self.wirkp.append(float(k_str_list[3]))
                    sum_weight += self.wirkp[-1]
                self.irkp = np.array(self.irkp)
                self.wirkp = np.divide(np.array(self.wirkp), sum_weight)
                continue

        # NBANDS, NGX/Y/Z
            if line.startswith('Dimension of arrays'):
                self.nbands = int(self.outlines[i+1].split()[-1])
                self.ngxyz  = [int(x) for x in self.outlines[i+6].split()[4:10:3]]
                self.ngfxyz  = [int(x) for x in self.outlines[i+7].split()[3:7:2]] 
                self.natoms = int(self.outlines[i+2].split()[-1])
        # Maximum number of planewaves
            if line.startswith('maximum number of plane-waves'):
                self.mplw = conv_string(line, int, -1)
        # ISPIN
            if line.startswith('ISPIN'):
                self.ispin = conv_string(line, int, 2)
        # NELECT: number of electrons
            if line.startswith('NELECT'):
                self.nelect = conv_string(line, int, 2)
        # ISMEAR and SIGMA: broadening information
            if line.startswith('ISMEAR'):
                self.ismear = conv_string(line.replace(';',''), int, 2)
                self.sigma  = conv_string(line.replace(';',''), float, 5)

        # lattice constants and inner coordinates 
        self.lattice, self.innerpos = self.get_pos(ionstep=0)


    def __check_finished(self):
        '''
        Check if the calculation has finished by checking whether the last line of outlines
        tells about the time summary
        '''
        if self.outlines[-1].strip().startswith('Voluntary context switches'):
            self.flag_finish = True
            return True
        else:
            self.flag_finish = False
            return False


    def __divide_ion_iteration(self):
        '''
        Check the line number of each ionic step
        '''
        ln_iteration = []
        i_ionstep = []
        self.flag_static = True

        for i in range(len(self.outlines)):
            words = self.outlines[i].split()
            if len(words) == 5:
                if words[1] == 'Iteration':
                    ln_iteration.append(i)
                    i_ionstep.append(int(words[2][:-1]))

        self.nisteps = i_ionstep[-1]
        if self.nisteps == 1:
            self.iterations = ln_iteration 
        else:
            self.flag_static = False
            self.iterations = [[] for i in range(self.nisteps)]
            for i in range(len(ln_iteration)):
                self.iterations[i_ionstep[i]-1].append(ln_iteration[i])


    def __get_ionic_data_region(self, ionstep, all_elesteps=False):
        '''
        Return the region of data in self.outlines for a particular ionic step.

        Parameters:
            ionstep: int
                The index of ionic step. Default -1 to return the information regarding the final results.
                0 for the initial data.
            all_elesteps: bool
                flag for controlling whether the whole range of electronic step or only the last step is included.
        '''
        
        if all_elesteps:
            end_ele = 0
        else:
            end_ele = -1

        if ionstep == 0:
            st_line = 0
            ed_line = self.iterations[0][0]
        elif ionstep == -1 or ionstep == self.nisteps:
            st_line = self.iterations[-1][end_ele]
            ed_line = len(self.outlines)-1
        elif ionstep < -1 or ionstep > self.nisteps:
            raise ValueError("The requested ionic step is too large. Maximum %d" % self.nisteps)
        else:
            st_line = self.iterations[ionstep-1][end_ele]
            ed_line = self.iterations[ionstep][0]

        return st_line, ed_line


    def get_pos(self, ionstep=-1):
        '''
        Get the initial lattice structure for a particular ionic step. Default return the last structure.

        Parameters:
            ionstep: int
                The index of ionic step. Default -1 to return the final structure.
                0 for the initial structure, i.e. POSCAR
        '''
        
        st_line, ed_line = self.__get_ionic_data_region(ionstep, False)
        innerpos = []

        for i in range(st_line, ed_line):
            #print(self.outlines[i])
            if self.outlines[i].strip().startswith('direct'):
                lattice = np.array([
                                [float(x) for x in self.outlines[i+1].split()[:3]],
                                [float(x) for x in self.outlines[i+2].split()[:3]],
                                [float(x) for x in self.outlines[i+3].split()[:3]],
                                   ])
            if ionstep == 0:
                if self.outlines[i].strip().startswith('position of ions in'):
                    for atom in range(self.natoms):
                        innerpos.append([float(x) for x in self.outlines[i+1+atom].split()])
                    break
            else:
                if self.outlines[i].strip().startswith('POSITION'):
                    for atom in range(self.natoms):
                        innerpos.append([float(x) for x in self.outlines[i+2+atom].split()[:3]])
                    break

        innerpos = np.array(innerpos)
        # convert cartisian to direct coordinates
        if ionstep != 0:
            innerpos = np.dot(innerpos, np.linalg.inv(lattice))
            # the inner coordinates are rounded, which may cause numerical error
            innerpos = np.round(innerpos,5)
        return lattice, innerpos


    def get_total_force(self, iatom=None, ionstep=-1):
        '''
        Get the total-force information for a particular ionic step.

        Parameters:
            iatom: int
                The index of atom to output the total force
            ionstep: int
                The index of ionic step. Default -1 to return the final step.
                ionstep=0 will give the same result as ionstep=1.

        Returns:
            total_force: numpy array
                The total force on all or one atom. 
                If iatom is set as a non-negative value, shape(3). Otherwise shape(natoms, 3)
        '''
       
        if ionstep==0:
            st_line, ed_line = self.__get_ionic_data_region(1, False)
        else:
            st_line, ed_line = self.__get_ionic_data_region(ionstep, False)

        self.total_force = []

        for i in range(st_line, ed_line):
            #print(self.outlines[i])
            if self.outlines[i].strip().startswith('POSITION'):
                for atom in range(self.natoms):
                    self.total_force.append([float(x) for x in self.outlines[i+2+atom].split()[3:]])
                break

        self.total_force = np.array(self.total_force)

        if iatom is None:
            return self.total_force
        else:
            try:
                assert int(iatom) < self.natoms
            except ValueError:
                raise ValueError('invalid index of atom')
            except AssertionError:
                raise ValueError('iatom out of range (%d)' % self.natoms)
            else:
                return self.total_force[iatom]


    def get_band_structure(self, ionstep=-1):
        '''
        Get the band structure for a particular ionic step.

        Parameters:
            ionstep: int
                The index of ionic step. -1 to return the band structure of final structure
                0 for that of the initial structure
        '''
        st_line, ed_line = self.__get_ionic_data_region(ionstep, False)
        # find the starting line for the band structure
        for i in range(st_line,ed_line):
            line = self.outlines[i].strip()
            if line.startswith('E-fermi'):
                ln_e_fermi = i
                break

        nbands = self.nbands
        ispin = self.ispin
        nirkp = self.nirkp

        kpt_block = nbands + 3
        spin_block = kpt_block * nirkp + 2

        self.efermi = conv_string(self.outlines[ln_e_fermi], float, 2)
        self.eigenene = []
        self.occ = []
        vb = self.nelect / 2 - 1
        cb = self.nelect / 2 
        self.vbm = -10000.0
        self.cbm =  10000.0
        self.vbm_loc =  [0, 0]
        self.cbm_loc =  [0, 0]
        self.Eg_dir = 10000.0
        self.dir_loc =  [0, 0]

        if ispin == 1:
            ln_spin_1_k_1 = ln_e_fermi + 3
        else:
            ln_spin_1_k_1 = ln_e_fermi + 5 # extra "spin component" and empty lines

        for spin in range(ispin):
            e_spin = []
            occ_spin = []
            for kpt in range(nirkp):
                e_kpt = []
                occ_kpt = []
                band_st_line = ln_spin_1_k_1 + spin*spin_block + kpt*kpt_block + 2
                for line in self.outlines[band_st_line:band_st_line+nbands]:
                    words = line.split()
                    e_kpt.append(float(words[1]))
                    occ_kpt.append(float(words[2]))

                Eg_dir = e_kpt[cb] - e_kpt[vb]
                if Eg_dir < self.Eg_dir:
                    self.Eg_dir  = Eg_dir
                    self.dir_loc = [spin, kpt]

                if e_kpt[vb] > self.vbm:
                    self.vbm     = e_kpt[vb]
                    self.vbm_loc = [spin, kpt]
                if e_kpt[cb] < self.cbm:
                    self.cbm     = e_kpt[cb]
                    self.cbm_loc = [spin, kpt]

                e_spin.append(e_kpt)
                occ_spin.append(e_kpt)
            self.eigenene.append(e_spin)
            self.occ.append(occ_spin)

        # the fundamental band gap
        self.Eg = self.cbm - self.vbm
        self.Eg_vbm = self.eigenene[self.vbm_loc[0]][self.vbm_loc[1]][cb] - self.eigenene[self.vbm_loc[0]][self.vbm_loc[1]][vb]
        self.Eg_cbm = self.eigenene[self.cbm_loc[0]][self.cbm_loc[1]][cb] - self.eigenene[self.cbm_loc[0]][self.cbm_loc[1]][vb]
        self.eigenene = np.array(self.eigenene)
        self.occ = np.array(self.occ)

        return self.eigenene, self.occ


    def get_gap(self, vb=None, cb=None):
        '''
        Print the band gap information. Also return the fundamental band gap.

        Parameters:
            vb: int
                the index of valence band (>=1)
            cb: int
                the index of conduction band (>=1)
        '''

        if vb is None and cb is None:
            return self.Eg
        elif vb is None:
            vb = self.nelect / 2
        elif cb is None:
            cb = self.nelect / 2 + 1

        # set the valence and 
        # if the band structure has not been extracted, use get_band_structure to obtain that of the final structure
        if self.eigenene is None:
            self.get_band_structure()


def read_fermi(outcar='OUTCAR'):

    fermi_level = 0.0

    if os.path.exists(outcar):
        try:
            fermi_level = float(sp.check_output("awk '/E-fermi/ {print $3}' %s" % outcar, shell=True))
            print("Fermi level found: %6.4f" % fermi_level)
        except ValueError:
            print("WARNING: error in reading E-fermi from %s. Fermi level is set to 0." % outcar)
    else:
        print("WARNING: %s is not found. Fermi level is set to 0." % outcar)

    return fermi_level


def vasp_anal_get_outcar(keyin,index=-1,outcar='OUTCAR'):
    '''
    return the value of key in outcar.
    '''
    key = keyin.lower()
# maximum number of plane wave, i.e. the maximum size of representation matrix
    if key=='mnpw':
        mnpw = sp.check_output("awk '/maximum number of/ {print $5}' %s | tail -1" % outcar,shell=True)
        return int(mnpw)
# NBANDS
    if key in ['nb', 'nbands']:
        nb = sp.check_output("awk '/NBANDS/ {print $15}' %s | head -1" % outcar,shell=True)
        return nb
# ENCUT
    if key in ['encut']:
        encut = sp.check_output("awk '/ENCUT/ {print $3}' %s | head -1" % outcar,shell=True)
        encut = int(float(encut))
        return encut
# converged G.S. energy
    if key in ['ene', 'energy']:
        ene = sp.check_output("awk '/without/ {print $7}' %s | tail -1" % outcar,shell=True)
        ene = float(ene)
        return ene
# total number of irreducible k-points
    if key in ['nkp']:
        nkp = sp.check_output("awk '/NKPTS/ {print $4}' %s | head -1" % outcar,shell=True)
        nkp = int(nkp)
        return nkp
    if key in ['efermi', 'e-fermi', 'fermi']:
        efermi = float(sp.check_output("awk '/E-fermi/ {print $3}' %s | tail -1" % outcar, shell=True))
        return efermi
