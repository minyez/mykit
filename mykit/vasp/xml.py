# coding = utf-8
'''
'''
class vasprunxmlError(Exception):
    pass


class vasprunxml:

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
    def sum_atom_l_comp(self, spin, band, kp, at_index, pw_index):
        weigh = 0
        for at in at_index:
            for pw in pw_index:
#                weigh += self.pwdata[spin][band][kp][at][pw]
                weigh += self.pwdata[spin][band][kp][at][pw]
        return weigh