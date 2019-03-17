# coding=utf-8

# ====================================================
#     File Name : pv_classes.py
# Creation Date : 30-10-2017
#    Created By : Min-Ye Zhang
#       Contact : stevezhang@pku.edu.cn
#       Purpose : provide vasp classes for analysis
# ====================================================

from __future__ import print_function

import os
import sys
from xml.etree import ElementTree as etree

import numpy as np

from mykit.core.utils import common_ss_conv

# ====================================================

class vasp_read_outcar:
    
    def __init__(self, OutFile='OUTCAR', verbose=True):
        
        self.filename = OutFile
        self.eigenene = None
        self.occ = None
        self.verbose = verbose

        # self.__print(" Reading OUTCAR: %s" % OutFile)
        with open(OutFile,'r') as f:
            self.outlines = f.readlines()
        self.__divide_ion_iteration()
        self.__init_calc_params()
        self.__check_finished()
        # if not self.flag_finish:
        #     common_print_warn("Task has not finished.", func_level=0, verbose=self.verbose)

    # def __print(self, pstr):
    #     common_print_verbose_bool(pstr, self.verbose)

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
                self.mplw = common_ss_conv(line, -1, int)
        # ISPIN
            if line.startswith('ISPIN'):
                self.ispin = common_ss_conv(line, 2, int)
        # NELECT: number of electrons
            if line.startswith('NELECT'):
                self.nelect = common_ss_conv(line, 2, int)
        # ISMEAR and SIGMA: broadening information
            if line.startswith('ISMEAR'):
                self.ismear = common_ss_conv(line.replace(';',''), 2, int)
                self.sigma  = common_ss_conv(line.replace(';',''), 5, float)

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

        self.efermi = common_ss_conv(self.outlines[ln_e_fermi], 2, float)
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

        

# ====================================================

# class vasp_read_poscar:
#     def check_mcenter(self):
#         '''
#         give the coordinate of the geometric center
#         '''
#         pass


#     def action_centering(self,zdirt1=3,zdirt2=None,zdirt3=None):
#         self.__print(" In action_centering:")
#         dirt_list = [zdirt1]
#         if (zdirt2 is not None) and (zdirt2 not in dirt_list):
#             dirt_list.append(zdirt2)
#         if (zdirt3 is not None) and (zdirt3 not in dirt_list):
#             dirt_list.append(zdirt2)
#         if self.coor_type == 'cart':
#             self.__print(" Cartisian coordinate system detected. Centering is not supported yet. Pass")
#             return

#         for zdirt in dirt_list:
#             iz = zdirt - 1
#             if self.check_vacuum_pos(zdirt):
#                 self.__print("  - Vacuum in the middle detected. Not supported currently. Pass.")
#                 continue
#             else:
#                 surf_atom = [self.check_extreme_atom(0.0,zdirt,False,1.0),self.check_extreme_atom(0.0,zdirt,True,1.0)]
#                 # debug
#                 # print surf_atom
#                 shift = 0.5 - sum([self.innerpos[i-1][iz] for i in surf_atom])/2.0
#                 self.action_shift(shift,zdirt)
#         self.__print(" Complete centering.")
#         #self.write_innerpos_to_lines()


# ===========================================================
# slab-related utilities
# ===========================================================

#     def check_extreme_atom(self,target=0.0,zdirt=3,l_far=True,terminal=1.0):
#         if self.coor_type == 'cart':
#             self.__print(" In check_extreme_atom:")
#             self.__print("  - Cartisian coordinate system detected. Currently not supported. Pass")
#             return None
#         for x in [target,terminal]:
#             if x < 0.0 or x > 1.0:
#                 self.__print(" In check_extreme_atom:")
#                 self.__print("  - Invalid target or terminal in extreme check. Pass.")
#                 return None
#         iz = zdirt - 1
#         ia_far = 0
#         ia_clo = 0
#         dmax = 0.0
#         dmin = 1.0
#         for ia in range(self.natoms):
#             dist = self.innerpos[ia][iz]-target
#             if ( dist *(self.innerpos[ia][iz]-terminal)>0):
#                 continue
#             dist = abs(dist)
#             if dist > dmax:
#                 dmax = dist
#                 ia_far = ia
#             if dist < dmin:
#                 dmin = dist
#                 ia_clo = ia

#         if l_far:
#             return ia_far+1
#         else:
#             return ia_clo+1


#     def check_vacuum_pos(self,zdirt):
#         thres_vac = 0.02
#         iz = zdirt - 1
#         zmax = max([self.innerpos[i][iz] for i in range(self.natoms)])
#         zmin = min([self.innerpos[i][iz] for i in range(self.natoms)])
#         # if the vacuum is in the middle of the slab/wire, return True
#         # else False
#         if zmin < thres_vac or (1.0-zmax) < thres_vac:
#             return True
#         else:
#             return False


#     def action_shift(self,shift,zdirt):
#         self.__print("  - overall shift %6.4f in direction %d" % (shift,zdirt))
#         for ia in range(self.natoms):
#             self.action_single_atom_shift(ia+1,shift,zdirt)


#     def action_single_atom_shift(self,iatom,shift,zdirt):
#         iz = zdirt - 1
#         ia = iatom - 1
#         self.innerpos[ia][iz] = self.innerpos[ia][iz] + shift
#         if self.coor_type == 'dirt':
#             self.check_pbc(iatom)


#     def action_add_vacuum(self,vacadd,zdirt=3):
#         '''
#         add length to the vacuum region of slab model. In Angstrom unit.
#         1D case and middle vacuum case to be implemented
#         '''
#         self.__print(" In action_add_vacuum:")
#         #if abs(vacadd) < 0.1:
#         #    print "  - too small change in vacuum (abs <= 1.0 A). Pass."
#         #    return
#         # check the slab model
#         iz = zdirt - 1
#         for ix in range(3):
#             if (iz is not ix) and self.lattice[iz][ix] > 1.0E-2:
#                 return
#             if (iz is not ix) and self.lattice[ix][iz] > 1.0E-2:
#                 return

#         if self.check_vacuum_pos(zdirt):
#             self.__print("  - Vacuum in the middle detected. Not supported currently. Pass.")
#             return

#         if self.coor_type == 'cart':
#             self.action_cart2dirt(False) # switch to direct system
#         zmax = max([self.innerpos[i][iz] for i in range(self.natoms)])
#         zmin = min([self.innerpos[i][iz] for i in range(self.natoms)])
#         vac_ori =  self.lattice[iz][iz] * (1.0 - (zmax-zmin))
#         self.__print("  - Original vacuum thickness: %8.5f" % vac_ori)

#         self.action_dirt2cart(False)     # switch to cartisian

#         self.lattice[iz][iz] = self.lattice[iz][iz] + vacadd
#         # shift the atoms in zdirt by half of vacadd
#         self.action_shift(vacadd/2.0,zdirt,False)
#         self.action_cart2dirt(True)

#         #self.action_centering(zdirt)
#         self.write_lattice_to_lines()


#     def check_sort_index(self,zdirt=3):
#         iz = zdirt - 1
#         Coords_all = [self.innerpos[i][iz] for i in range(self.natoms)]
#         sorted_index = sorted(range(self.natoms),key=lambda k:Coords_all[k])
#         return sorted_index


#     def action_sort_coor(self,zdirt=3):
#         '''
#         Sort atoms of each type
#         '''
#         # skip the heading, save the starting line of each atom block
#         iz = zdirt - 1
#         blocka = [0]
#         for i in range(self.ntypes-1):
#             blocka.append(blocka[i]+self.atom_num[i])

#         for i in range(self.ntypes):
#             Coords = [self.innerpos[x][iz] for x in range(blocka[i],blocka[i]+self.atom_num[i])]
# #            print Coords
#             index = sorted(range(self.atom_num[i]),key=lambda k:Coords[k])
# #            print index
#             temp_list = []
#         # sort the coordinates for each type of atoms
#             for j in range(self.atom_num[i]):
#                 temp_list.append(self.innerpos[blocka[i]+index[j]])
#             self.innerpos[blocka[i]:blocka[i]+self.atom_num[i]] = temp_list

#         self.write_innerpos_to_lines()
#         return self.check_sort_index(zdirt)

# ====================================================

class vasp_read_xml():

    '''
    Class to read and analyse the data from the vasprun.xml
    '''

    def __init__(self,agrv=[], verbose=True):
        self.verbose = verbose
        tree = etree.parse('vasprun.xml')
        self.root = tree.getroot()
        self.init_section()
        self.read_atominfo()
        self.read_para()
        self.read_klist()
        self.read_eigen()
        try:
            if 'pw' in agrv:
                self.read_pwdata()
        except AttributeError:
            self.__print("Warning: no PDOS data found")

    def __print(self, pstr):
        # common_print_verbose_bool(pstr, self.verbose)
        print(pstr, self.verbose)


    def init_section(self):
        self.para = self.root.find('parameters')
        # the last calculation which has the final eigenvalues
        self.calc = self.root.findall('calculation')[-1]
        self.atominfo = self.root.find('atominfo')
        self.kps = self.root.find('kpoints')


    def read_para(self):
        # ISPIN: root->parameter->separator name='electronic'->separator name='elecronic spin'->[0]
        # NBANDS: root->parameter->separator name='electronic'->i name='NBANDS'
        para_e = self.para.find('.//separator[@name="electronic"]')
        spin = para_e.find('.//separator[@name="electronic spin"]')
        self.nelec = int(float(para_e.find('.//i[@name="NELECT"]').text))
        self.ispin = int(spin[0].text)
        self.nbands = int(para_e.find('.//i[@name="NBANDS"]').text)


    def read_pwdata(self):
        # partial wave strings: root->dos->partial->array->field[1:] # the first is energy label
        partial_array = self.calc.find('dos').find('partial').find('array')
        self.str_pwaves = [pwave.text.strip() for pwave in partial_array.findall('field')[1:]]
        # pdos data: root->calculation->projected->array
        proj =  self.calc.find('projected')
        pwav_dataset = proj.findall('array')[0][-1] # contains ISPIN dataset(s)

        # initialize the pwav data
        self.pwdata = np.zeros([self.ispin,self.nkp,self.nbands,self.natom,len(self.str_pwaves)])

        for spin in range(self.ispin):
            for kp in range(self.nkp):
                for band in range(self.nbands):
                    for atom in range(self.natom):
                        data = np.array([float(x) for x in pwav_dataset[spin][kp][band][atom].text.split()])


    def read_atominfo(self):
        self.natom = int(self.atominfo[0].text)
        self.ntype = int(self.atominfo[1].text)
        self.atoms = [0]*self.ntype
        self.type = []
        self.type_list = []
        type_list = []
        atomdata = self.atominfo.find('array').find('set')
        for atom in atomdata:
            type_name = atom[0].text.split()[0]
            self.type_list.append(type_name)
            if type_name not in self.type:
                self.type.append(type_name)
            self.atoms[int(atom[1].text)-1] += 1


    def read_eigen(self):
        self.eigen = np.zeros([self.ispin,self.nkp,self.nbands])
        eigen_data = self.calc.find('eigenvalues').find('array').find('set')
        for spin in range(self.ispin):
            for kp in range(self.nkp):
                self.eigen[spin,kp] = np.array([float(x.text.split()[0]) for x in eigen_data[spin][kp]])


    def get_gap(self):
        # return the gap from eigenvalue information
        ecbm = 100000.0
        evbm = -100000.0
        vbm = self.nelec/2
        cbm = self.nelec/2 + 1
        for spin in range(self.ispin):
            for kp in range(self.nkp):
#                print self.eigen[spin,kp,vbm-1],self.eigen[spin,kp,cbm-1]
                if (self.eigen[spin,kp,vbm-1]> evbm):
                    evbm = self.eigen[spin,kp,vbm-1]
                if (self.eigen[spin,kp,cbm-1]< ecbm):
                    ecbm = self.eigen[spin,kp,cbm-1]
        gap = ecbm - evbm
        return gap


    def read_klist(self):
#       klist_index is 1 if auto generator is used
#       or 0 if mannually included
        if self.kps.find('generation') is not None:
            ki = 1
        else:
            ki = 0
        self.kplist = [[ float(kvec) for kvec in kp.text.split()] for kp in self.kps[ki]]
        self.kpweigh = np.array([float(kp.text) for kp in self.kps[ki+1]])
        self.nkp = len(self.kplist)


    def atoms_index(self,atom_type):
        try:
            itype = int(atom_type)
            assert itype <= self.ntype
        except:
            return None
        if atom_type == 0:
            return list(range(self.natom))
        elif atom_type == 1:
            return list(range(self.atoms[0]))
        else:
            list1 = list(range(sum(self.atoms[:itype-1])))
            list2 = list(range(sum(self.atoms[:itype])))
            return [x for x in list2 if x not in list1]


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
