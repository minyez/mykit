#!/usr/bin/env python
# coding=utf-8

# ====================================================
#     File Name : pv_classes.py
# Creation Date : 30-10-2017
#    Created By : Min-Ye Zhang
#       Contact : stevezhang@pku.edu.cn
#       Purpose : provide vasp classes for analysis
# ====================================================

import sys
import numpy as np
from xml.etree import ElementTree as etree

class vasp_read_poscar:
    def __init__(self,PosFile='POSCAR'):
        print " Reading POSCAR: %s" % PosFile
        with open(PosFile,'r') as f:
            self.poslines = f.readlines()
        self.__comment_info()
        self.__latt_info()
        self.__atom_info()
        self.__coor_info()
        print " POSCAR read."


    def __comment_info(self):
        self.comment = self.poslines[0]   # the comment line


    def __latt_info(self):
        # unscaled lattice vectors
        # scaled lattice volume
        print " Loading lattice info..."
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
        print " Loading atomic info..."
        posls = self.poslines
        self.atom_type = posls[5].split()
        self.ntypes = len(self.atom_type)
        self.atom_num = [int(i) for i in posls[6].split()]
        self.natoms = sum(self.atom_num)
        self.composition = ""
        for i in xrange(len(self.atom_type)):
            self.composition += self.atom_type[i]+str(self.atom_num[i])


    def __coor_info(self):
        print " Loading inner coordinates..."
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
            print " Unknown coordinate system tag. Exit"
            sys.exit(1)
        self.poslines[7] = 'Selective dynamics\n'+coor_type
        self.write_innerpos_from_lines()
        print " Checking periodic boundary condition..."
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
            print " In check_pbc:"
            print "  - Cartisian system detected. PBC check not supported yet. Pass"
            return
        else:
            for iz in xrange(3):
                while (self.innerpos[ia][iz] >= 1.0):
                    print "  - PBC check. Atom %3d, Coord. %d shift down" % (iatom,iz+1)
                    self.innerpos[ia][iz] = self.innerpos[ia][iz] - 1.0
                while (self.innerpos[ia][iz] < 0.0):
                    print "  - PBC check. Atom %3d, Coord. %d shift up" % (iatom,iz+1)
                    self.innerpos[ia][iz] = self.innerpos[ia][iz] + 1.0


    def check_atomtype(self,iatom):
        if iatom <= 0 or iatom > self.natoms:
            print "Invalid atom index (should be positive and <= natoms)"
            return None
        ia = iatom
        for itype in xrange(self.ntypes):
            if ia <= self.atom_num[itype]:
                return self.atom_type[itype]
            else:
                ia = ia - self.atom_num[itype]


    def check_mcenter(self):
        '''
        give the coordinate of the center of mass
        '''
        pass

    def action_centering(self,zdirt1=3,zdirt2=None,zdirt3=None):
        dirt_list = [zdirt1-1]
        if (zdirt2 is not None) and ((zdirt2-1) not in dirt_list):
            dirt_list.append(zdirt2-1)
        if (zdirt3 is not None) and ((zdirt3-1) not in dirt_list):
            dirt_list.append(zdirt2-1)
        if self.coor_type == 'cart':
            print " Cartisian coordinate system detected. Centering is not supported yet. Pass"
            return

        for iz in dirt_list:
            if self.check_vacuum_pos(iz+1):
                print "  - Vacuum in the middle detected. Not supported currently. Pass."
                continue
            else:
                surf_atom = [self.check_extreme_atom(0.0,3,False,1.0),self.check_extreme_atom(0.0,3,True,1.0)]
                shift = 0.5 - sum([self.innerpos[i-1][iz] for i in surf_atom])/2.0
                if shift >= 1E-3:
                    self.action_shift(shift,iz+1)
        self.write_innerpos_to_lines()


    def action_dirt2cart(self,l_write=True):
        if self.coor_type == 'cart':
            print " Cartisian coordinate system detected. Nothing to do."
            return
        for ia in xrange(self.natoms):
            self.innerpos[ia] = np.dot(self.lattice,self.innerpos[ia])
        self.coor_type = 'cart'
        if l_write:
            self.write_coortype_to_lines()
            self.write_innerpos_to_lines()


    def action_cart2dirt(self,l_write=True):
        if self.coor_type == 'dirt':
            print " Direct coordinate system detected. Nothing to do."
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


    def write_poscar(self,OutFile='POSCAR_new',l_reload=False):
        print " Write POSCAR: %s" % OutFile
        # reload all objects to poslines
        if l_reload:
            self.action_change_scale(self.scale)
            self.write_coortype_to_lines()
            self.write_innerpos_to_lines()
            self.write_lattice_to_lines()
        with open(OutFile,'w+') as f:
            for i in xrange(len(self.poslines)):
                f.write(self.poslines[i])


# ===========================================================
# slab utilities
# ===========================================================

    def check_extreme_atom(self,target=0.0,zdirt=3,l_far=True,terminal=1.0):
        if self.coor_type == 'cart':
            print " In check_extreme_atom:"
            print "  - Cartisian coordinate system detected. Currently not supported. Pass"
            return None
        for x in [target,terminal]:
            if x < 0.0 or x > 1.0:
                print " In check_extreme_atom:"
                print "  - Invalid target or terminal in extreme check. Pass."
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


    def action_shift(self,shift,zdirt,verbose=True):
        if verbose:
            print "  - overall shift %6.4f in direction %d" % (shift,zdirt)
        for ia in xrange(self.natoms):
            self.action_single_atom_shift(ia+1,shift,zdirt)


    def action_single_atom_shift(self,iatom,shift,zdirt):
        iz = zdirt - 1
        ia = iatom - 1
        self.innerpos[ia][iz] = self.innerpos[ia][iz] + shift
        if self.coor_type == 'dirt':
            self.check_pbc(iatom)


    def action_add_vacuum(self,vacadd,zdirt=3):
        # add length to the vacuum region of slab model. In Angstrom unit
        # 1D case and middle vacuum case to be implemented
        if abs(vacadd) < 1.0:
            print " In action_add_vacuum:"
            print "  - too small change in vacuum (abs <= 1.0 A). Pass."
            return
        # check the slab model
        iz = zdirt - 1
        for ix in xrange(3):
            if (iz is not ix) and self.lattice[iz][ix] > 1.0E-2:
#                print
                return
            if (iz is not ix) and self.lattice[ix][iz] > 1.0E-2:
                return

        if self.check_vacuum_pos(zdirt):
            print " In action_add_vacuum:"
            print "  - Vacuum in the middle detected. Not supported currently. Pass."
            return

        if self.coor_type == 'cart':
            self.action_cart2dirt(False) # switch to direct system
        zmax = max([self.innerpos[i][iz] for i in xrange(self.natoms)])
        zmin = min([self.innerpos[i][iz] for i in xrange(self.natoms)])
        vac_ori =  self.lattice[iz][iz] * (1.0 - (zmax-zmin))
        print " In action_add_vacuum:"
        print "  - Original vacuum thickness: %8.5f" % vac_ori
        self.action_centering(zdirt)
        self.action_dirt2cart(False)     # switch to cartisian

        self.lattice[iz][iz] = self.lattice[iz][iz] + vacadd
        self.action_shift(vacadd/2.0,zdirt)
        self.action_cart2dirt(True)
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


# ===========================================================

class vasp_read_xml():

    def __init__(self,agrv=[]):
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
            print "Warning: no PDOS data found"

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



