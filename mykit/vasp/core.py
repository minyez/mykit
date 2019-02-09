#!/usr/bin/env python
# coding=utf-8

# ====================================================
#     File Name : pv_classes.py
# Creation Date : 30-10-2017
#    Created By : Min-Ye Zhang
#       Contact : stevezhang@pku.edu.cn
#       Purpose : provide vasp classes for analysis
# ====================================================

from __future__ import print_function
import sys,os
import numpy as np
from xml.etree import ElementTree as etree
from mykit.common.utils import print_warn, ss_conv, print_verbose_bool

class vasp:

class incar:

    def __init__(self, InFile=None, \
                 encut=None, xc=None, lreal='.FALSE.', \
                 ispin=1, istart=0, ediff=1E-7, \
                 verbose=True):

        self.__init_tags_dict() 
        self.readinfile = InFile
        self.verbose = verbose
        self.header = None
        self.xc = xc

        self.incar_tags_dict = {}
        if InFile is not None and os.path.exists(InFile):
            self.__read_existing_tags(InFile)
        self.mod_xc_tags(xc)
        self.add_single_tag('ISTART', istart)
        self.add_single_tag('ISPIN', ispin)
        self.add_single_tag('LREAL', lreal)
        self.add_single_tag('ENCUT', encut)

    def __print(self, pstr):
        common_print_verbose_bool(pstr, self.verbose)

    def __getitem__(self, key):
        return self.incar_tags_dict[key]


    def __init_tags_dict(self):
        self.bool_tags  = ['LCHARG','LWAVE','LDIPOL','LVHAR','LVTOT','LHFCALC','LASPH','LMIXTAU']
        self.int_tags   = ['NGX','NGY','NGZ','NGXF','NGYF','NGZF','NBANDS','ISTART',\
                           'ICHARG','ISPIN','NWRITE','NELM','NSW','IBRION','ISIF','ISYM',\
                           'IALGO','ISMEAR','NKRED','NKREDX','NKREDY','NKREDZ','ODDONLY',\
                           'EVENONLY','LMAXTAU','IDIPOL','NELMIN','NPAR']
        self.float_tags = ['AEXX','AGGAX','AGGAC','ALDAC','ENCUT','ENCUTGW','EDIFF',\
                           'EDIFFG','POTIM','NELECT','SMASS','AMIX','BMIX','SIGMA',\
                           'HFSCREEN','TIME']
        self.str_tags   = ['PREC','PRECFOCK','ALGO','LREAL','GGA','METAGGA']
        self.long_tags  = ['MAGMOM','DIPOL']

        self.all_tags = self.bool_tags + self.int_tags + \
                        self.float_tags + self.str_tags + self.long_tags

        self.tags_conv_dict = {}
        for tag in self.bool_tags: 
            self.tags_conv_dict.update({tag: bool})
        for tag in self.int_tags: 
            self.tags_conv_dict.update({tag: int})
        for tag in self.float_tags: 
            self.tags_conv_dict.update({tag: float})
        for tag in self.long_tags:
            self.tags_conv_dict.update({tag: str})
        for tag in self.str_tags:
            self.tags_conv_dict.update({tag: str})

        self.tag_PBE0 = {'GGA':'PE', 'ALGO':'ALL', 'LHFCALC':'.TRUE.', \
                    'TIME':0.4, 'PRECFOCK':'NORMAL'}
        self.tag_HSE06 = self.tag_PBE0.copy()
        self.tag_HSE06.update({'HFSCREEN':0.2})
        self.tag_HF = self.tag_PBE0.copy() 
        self.tag_HF.update({'AEXX':1.0})

        self.xc_dict = {
                    'LDA'    : {'GGA'     : 'CA'} ,
                    'PBE'    : {'GGA'     : 'PE'} ,
                    'PBESOL' : {'GGA'     : 'PS'} ,
                    'RPBE'   : {'GGA'     : 'RP'} ,
                    'SCAN'   : {'METAGGA' : 'SCAN', 'LASPH' : '.TRUE.'} ,
                    'PBE0'   : self.tag_PBE0  ,
                    'HSE06'  : self.tag_HSE06 ,
                    'HF'     : self.tag_HF
                  }


    def __read_existing_tags(self, InFile):

        with open(InFile, 'rw') as h_InFile:
            self.inlines = h_InFile.readlines()

        # check if the INCAR file is empty
        if len(self.inlines) == 0:
            return
        else:
            self.__print(" Read input tags from %s" % InFile)

        # if = is not found in the first line, set it as the INCAR header
        # or it is just a tag string
        if '=' not in self.inlines[0]:
            self.header = self.inlines[0].strip()
            self.inlines = self.inlines[1:]
        
        for iline in self.inlines:
            line_striped = iline.strip()
            # skip for comment and empty line
            if line_striped.startswith("#") or len(line_striped) == 0:
                continue

            # first separate with colon for multiple tags
            tag_sep_list = line_striped.split(';')
            
            for tag in tag_sep_list:
                tag_words = [x.strip() for x in tag.split('=')]
                # check if comment is found in the value part.
                # if so, delete the comment
                if '#' in tag_words[1]:
                    tag_words[1] = tag_words[1].split('#')[0].strip()
                if tag_words[0] in self.tags_conv_dict:
                    tag_words[1] = self.tags_conv_dict[tag_words[0]](tag_words[1])
                self.incar_tags_dict.update({tag_words[0]:tag_words[1]})


    def add_single_tag(self, tag, value):
        if value in [None,'']:
            return
        self.incar_tags_dict.update({tag.upper(): value})


    def write_incar_file(self, InCARout='INCAR', replace=True):
        if os.path.exists(InCARout):
            # backup the old INCAR
            if replace:
                os.rename(InCARout, InCARout+'.old')
        with open(InCARout,'w') as h_incaro:
            if self.header is not None:
                h_incaro.write(self.header+'\n')
            # print the key in a sorted order
            for tag in sorted(self.incar_tags_dict.keys()):
                h_incaro.write(' %-8s = %s\n' % (tag, str(self.incar_tags_dict[tag])))


    def mod_xc_tags(self, xc):
        '''
        Modify XC related tags to the tag dictionary.
        '''
        if xc is None:
            pass
        else:
            if xc.upper() == self.xc:
                self.__print(" XC is already %s " % self.xc)
            else:
                if self.xc is not None:
                    for tag in self.xc_dict[self.xc].keys():
                        self.incar_tags_dict.pop(tag)
                
                xc_origin = self.xc
                self.xc = xc.upper()

                self.__print(" XC set to %6s, from %6s" % (self.xc, xc_origin))

                for tag in self.xc_dict[self.xc].keys():
                    self.incar_tags_dict.update({tag: self.xc_dict[self.xc][tag]})


        return

# ====================================================

class vasp_read_wavecar:

    def __init__(self, WavFile='WAVECAR', verbose=True):
        
        self.filename = WavFile
        self.verbose = verbose
        with open(WavFile,'rb') as h_WavFile:
            #self.wavlines = h_WavFile.readlines()
            pass

# ====================================================

class vasp_read_outcar:
    
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
            common_print_warn("Task has not finished.", func_level=0, verbose=self.verbose)

    def __print(self, pstr):
        common_print_verbose_bool(pstr, self.verbose)

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

class vasp_read_poscar:
    def __init__(self,PosFile='POSCAR', verbose=True):

        print(" Reading POSCAR: %s" % PosFile)
        with open(PosFile,'r') as f:
            self.poslines = f.readlines()
        self.filename = PosFile
        self.verbose = verbose
        self.__comment_info()
        self.__latt_info()
        self.__atom_info()
        self.__coor_info()
        print(" POSCAR read.")

    def __print(self, pstr):
        common_print_verbose_bool(pstr, self.verbose)

    def __save_to_lines(self):
        '''
        save poscar data to poscar lines
        '''
        # reload all objects to poslines
        self.action_change_scale(self.scale)
        self.write_coortype_to_lines()
        self.write_innerpos_to_lines()
        self.write_lattice_to_lines()


    def __comment_info(self):
        self.comment = self.poslines[0]   # the comment line


    def __latt_info(self):
        # unscaled lattice vectors
        # scaled lattice volume
        self.__print(" Loading lattice info...")
        posls = self.poslines
        self.scale = float(posls[1].split()[0])
        self.lattice = []
        for i in [2,3,4]:
            self.lattice.append(np.array([float(x) for x in posls[i].split()]))
        self.lenlat = [np.sqrt(np.dot(x,x)) for x in self.lattice]
        self.lattice = np.array(self.lattice)
        self.lenlat = np.array(self.lenlat)
        self.volume = np.power(self.scale,3) * np.dot(np.cross(self.lattice[0],self.lattice[1]),self.lattice[2])


    def __atom_info(self):
        self.__print(" Loading atomic info...")
        posls = self.poslines
        self.atom_type = posls[5].split()
        self.ntypes = len(self.atom_type)
        self.atom_num = [int(i) for i in posls[6].split()]
        self.natoms = sum(self.atom_num)
        self.composition = ""
        for i in xrange(len(self.atom_type)):
            self.composition += self.atom_type[i]+str(self.atom_num[i])


    def __coor_info(self):
        self.__print(" Loading inner coordinates...")
        posls = self.poslines
        coor_type = ''
        # delete the 'selevtive dynamics' line
        if posls[7].startswith('S') or posls[7].startswith('s'):
            del self.poslines[7]
        if posls[7].startswith('C') or posls[7].startswith('c') \
            or posls[7].startswith('K') or posls[7].startswith('k'):
            self.coor_type = 'cart'
            coor_type = 'Cartisian\n'
        elif posls[7].startswith('D') or posls[7].startswith('d'):
            self.coor_type = 'dirt'
            coor_type = 'Direct\n'
        else:
            self.__print(" Unknown coordinate system tag. Exit")
            sys.exit(1)
        self.poslines[7] = 'Selective dynamics\n'+coor_type
        self.write_innerpos_from_lines()
        self.__print(" Checking periodic boundary condition...")
        for i in xrange(self.natoms):
            self.check_pbc(i+1)


    def check_pbc(self,iatom):
        '''
        Check if the innerpos of [iatom]th atom is within the unit cell
        If not, shift the innerpos according to the periodic boundary condition
        Currently ONLY support DIRECT coordinate system
        '''
        # ia is the index of [iatom]th atom in self.innerpos[:]
        ia = iatom - 1
        if self.coor_type == 'cart':
            self.__print(" In check_pbc:")
            self.__print("  - Cartisian system detected. PBC check not supported yet. Pass")
            return
        else:
            for iz in xrange(3):
                while (self.innerpos[ia][iz] >= 1.0):
                    self.__print("  - PBC check. Atom %3d, Coord. %d shift down" % (iatom,iz+1))
                    self.innerpos[ia][iz] = self.innerpos[ia][iz] - 1.0
                while (self.innerpos[ia][iz] < 0.0):
                    self.__print("  - PBC check. Atom %3d, Coord. %d shift up" % (iatom,iz+1))
                    self.innerpos[ia][iz] = self.innerpos[ia][iz] + 1.0


    def check_atomtype(self,iatom):
        if iatom <= 0 or iatom > self.natoms:
            self.__print("Invalid atom index (should be positive and <= natoms)")
            return None
        ia = iatom
        for itype in xrange(self.ntypes):
            if ia <= self.atom_num[itype]:
                return self.atom_type[itype]
            else:
                ia = ia - self.atom_num[itype]


    def check_mcenter(self):
        '''
        give the coordinate of the geometric center
        '''
        pass


    def action_centering(self,zdirt1=3,zdirt2=None,zdirt3=None):
        self.__print(" In action_centering:")
        dirt_list = [zdirt1]
        if (zdirt2 is not None) and (zdirt2 not in dirt_list):
            dirt_list.append(zdirt2)
        if (zdirt3 is not None) and (zdirt3 not in dirt_list):
            dirt_list.append(zdirt2)
        if self.coor_type == 'cart':
            self.__print(" Cartisian coordinate system detected. Centering is not supported yet. Pass")
            return

        for zdirt in dirt_list:
            iz = zdirt - 1
            if self.check_vacuum_pos(zdirt):
                self.__print("  - Vacuum in the middle detected. Not supported currently. Pass.")
                continue
            else:
                surf_atom = [self.check_extreme_atom(0.0,zdirt,False,1.0),self.check_extreme_atom(0.0,zdirt,True,1.0)]
                # debug
                # print surf_atom
                shift = 0.5 - sum([self.innerpos[i-1][iz] for i in surf_atom])/2.0
                self.action_shift(shift,zdirt)
        self.__print(" Complete centering.")
        #self.write_innerpos_to_lines()


    def action_dirt2cart(self,l_write=True):
        if self.coor_type == 'cart':
            self.__print(" Cartisian coordinate system detected. Nothing to do.")
            return
        for ia in xrange(self.natoms):
            self.innerpos[ia] = np.dot(self.lattice,self.innerpos[ia])
        self.coor_type = 'cart'
        if l_write:
            self.write_coortype_to_lines()
            self.write_innerpos_to_lines()


    def action_cart2dirt(self,l_write=True):
        if self.coor_type == 'dirt':
            self.__print(" Direct coordinate system detected. Nothing to do.")
            return
        # use coordinate transformation x_cart = Ax_dirt, x_dirt = A^{-1}x_cart
        for ia in xrange(self.natoms):
            self.innerpos[ia] = np.dot(np.linalg.inv(self.lattice),self.innerpos[ia])
        self.coor_type = 'dirt'
        if l_write:
            self.write_coortype_to_lines()
            self.write_innerpos_to_lines()


    def action_change_scale(self,scale=1.0,l_write=True):
        self.scale = scale
        if l_write:
            self.poslines[1] = "%12.8f\n" % self.scale


    def write_lattice_to_lines(self):
        for i in xrange(3):
            self.poslines[i+2] = "%20.12f%20.12f%20.12f\n" % \
                           (self.lattice[i][0],self.lattice[i][1],self.lattice[i][2])


    def write_innerpos_from_lines(self):
        posls = self.poslines
        self.innerpos = []
        for n_line in xrange(8,8+self.natoms):
            self.innerpos.append(np.array([float(x) for x in posls[n_line].split()[0:3]]))


    def write_innerpos_to_lines(self):
        for ia in xrange(self.natoms):
            self.poslines[ia+8] =  "%17.10f%17.10f%17.10f" % \
                    (self.innerpos[ia][0],self.innerpos[ia][1],self.innerpos[ia][2]) \
                    + "  T  T  T " + '#%2s\n' % self.check_atomtype(ia+1)


    def write_coortype_to_lines(self):
        if self.coor_type == 'cart':
            self.poslines[7] = 'Selective dynamics\nCartisian\n'
        if self.coor_type == 'dirt':
            self.poslines[7] = 'Selective dynamics\nDirect\n'


    def write_poscar(self,OutFile='POSCAR_new',f_reload=True):
        self.__print(" Write POSCAR: %s" % OutFile)
        if f_reload:
            self.__print(" - Load modified atomic information to poslines")
            self.__save_to_lines()
        if OutFile is not None:
            if OutFile==self.filename:
                common_print_warn("OutFile same as input. Change input to %s.old" % self.filename, func_level=0, verbose=self.verbose)
                os.rename(self.filename,self.filename+'.old') 
                self.filename = self.filename + '.old'
            with open(OutFile,'w+') as f:
                for i in xrange(len(self.poslines)):
                    f.write(self.poslines[i])
            self.__print(" - Done output: %s" % OutFile)
        else:
            common_print_warn("Invalid output filename: None. Pass", func_level=0, verbose=self.verbose)

# ===========================================================
# slab-related utilities
# ===========================================================

    def check_extreme_atom(self,target=0.0,zdirt=3,l_far=True,terminal=1.0):
        if self.coor_type == 'cart':
            self.__print(" In check_extreme_atom:")
            self.__print("  - Cartisian coordinate system detected. Currently not supported. Pass")
            return None
        for x in [target,terminal]:
            if x < 0.0 or x > 1.0:
                self.__print(" In check_extreme_atom:")
                self.__print("  - Invalid target or terminal in extreme check. Pass.")
                return None
        iz = zdirt - 1
        ia_far = 0
        ia_clo = 0
        dmax = 0.0
        dmin = 1.0
        for ia in xrange(self.natoms):
            dist = self.innerpos[ia][iz]-target
            if ( dist *(self.innerpos[ia][iz]-terminal)>0):
                continue
            dist = abs(dist)
            if dist > dmax:
                dmax = dist
                ia_far = ia
            if dist < dmin:
                dmin = dist
                ia_clo = ia

        if l_far:
            return ia_far+1
        else:
            return ia_clo+1


    def check_vacuum_pos(self,zdirt):
        thres_vac = 0.02
        iz = zdirt - 1
        zmax = max([self.innerpos[i][iz] for i in xrange(self.natoms)])
        zmin = min([self.innerpos[i][iz] for i in xrange(self.natoms)])
        # if the vacuum is in the middle of the slab/wire, return True
        # else False
        if zmin < thres_vac or (1.0-zmax) < thres_vac:
            return True
        else:
            return False


    def action_shift(self,shift,zdirt):
        self.__print("  - overall shift %6.4f in direction %d" % (shift,zdirt))
        for ia in xrange(self.natoms):
            self.action_single_atom_shift(ia+1,shift,zdirt)


    def action_single_atom_shift(self,iatom,shift,zdirt):
        iz = zdirt - 1
        ia = iatom - 1
        self.innerpos[ia][iz] = self.innerpos[ia][iz] + shift
        if self.coor_type == 'dirt':
            self.check_pbc(iatom)


    def action_add_vacuum(self,vacadd,zdirt=3):
        '''
        add length to the vacuum region of slab model. In Angstrom unit.
        1D case and middle vacuum case to be implemented
        '''
        self.__print(" In action_add_vacuum:")
        #if abs(vacadd) < 0.1:
        #    print "  - too small change in vacuum (abs <= 1.0 A). Pass."
        #    return
        # check the slab model
        iz = zdirt - 1
        for ix in xrange(3):
            if (iz is not ix) and self.lattice[iz][ix] > 1.0E-2:
                return
            if (iz is not ix) and self.lattice[ix][iz] > 1.0E-2:
                return

        if self.check_vacuum_pos(zdirt):
            self.__print("  - Vacuum in the middle detected. Not supported currently. Pass.")
            return

        if self.coor_type == 'cart':
            self.action_cart2dirt(False) # switch to direct system
        zmax = max([self.innerpos[i][iz] for i in xrange(self.natoms)])
        zmin = min([self.innerpos[i][iz] for i in xrange(self.natoms)])
        vac_ori =  self.lattice[iz][iz] * (1.0 - (zmax-zmin))
        self.__print("  - Original vacuum thickness: %8.5f" % vac_ori)

        self.action_dirt2cart(False)     # switch to cartisian

        self.lattice[iz][iz] = self.lattice[iz][iz] + vacadd
        # shift the atoms in zdirt by half of vacadd
        self.action_shift(vacadd/2.0,zdirt,False)
        self.action_cart2dirt(True)

        #self.action_centering(zdirt)
        self.write_lattice_to_lines()


    def check_sort_index(self,zdirt=3):
        iz = zdirt - 1
        Coords_all = [self.innerpos[i][iz] for i in xrange(self.natoms)]
        sorted_index = sorted(range(self.natoms),key=lambda k:Coords_all[k])
        return sorted_index


    def action_sort_coor(self,zdirt=3):
        '''
        Sort atoms of each type
        '''
        # skip the heading, save the starting line of each atom block
        iz = zdirt - 1
        blocka = [0]
        for i in xrange(self.ntypes-1):
            blocka.append(blocka[i]+self.atom_num[i])

        for i in xrange(self.ntypes):
            Coords = [self.innerpos[x][iz] for x in xrange(blocka[i],blocka[i]+self.atom_num[i])]
#            print Coords
            index = sorted(range(self.atom_num[i]),key=lambda k:Coords[k])
#            print index
            temp_list = []
        # sort the coordinates for each type of atoms
            for j in xrange(self.atom_num[i]):
                temp_list.append(self.innerpos[blocka[i]+index[j]])
            self.innerpos[blocka[i]:blocka[i]+self.atom_num[i]] = temp_list

        self.write_innerpos_to_lines()
        return self.check_sort_index(zdirt)

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
        common_print_verbose_bool(pstr, self.verbose)


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

        for spin in xrange(self.ispin):
            for kp in xrange(self.nkp):
                for band in xrange(self.nbands):
                    for atom in xrange(self.natom):
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
        for spin in xrange(self.ispin):
            for kp in xrange(self.nkp):
                self.eigen[spin,kp] = np.array([float(x.text.split()[0]) for x in eigen_data[spin][kp]])


    def get_gap(self):
        # return the gap from eigenvalue information
        ecbm = 100000.0
        evbm = -100000.0
        vbm = self.nelec/2
        cbm = self.nelec/2 + 1
        for spin in xrange(self.ispin):
            for kp in xrange(self.nkp):
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
            return list(xrange(self.natom))
        elif atom_type == 1:
            return list(xrange(self.atoms[0]))
        else:
            list1 = list(xrange(sum(self.atoms[:itype-1])))
            list2 = list(xrange(sum(self.atoms[:itype])))
            return [x for x in list2 if x not in list1]


    def pwave_index(self,lcomponent):
        index_l = []
        # total wave
        if lcomponent == 't':
            index_l = list(xrange(len(self.str_pwaves)))
        for x in self.str_pwaves:
            if x.startswith(lcomponent):
                index_l.append(self.str_pwaves.index(x))
        return index_l


# sum the component of the atoms in at_index and partial waves in pw_index
# can be more pythonic
    def sum_atom_l_comp(self,spin,band,kp,at_index,pw_index):
        weigh = 0
        for at in at_index:
            for pw in pw_index:
#                weigh += self.pwdata[spin][band][kp][at][pw]
                weigh += self.pwdata[spin][band][kp][at][pw]
        return weigh


