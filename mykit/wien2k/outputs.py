# coding = utf-8
"""Uitlities for analysing WIEN2k output files
"""

import os
from collections import OrderedDict
from mykit.wien2k.utils import get_casename, read_lm_data
from mykit.core.utils import conv_string, find_str_matched
from mykit.core.numeric import Prec
from mykit.core.log import Verbose
from mykit.core.bandstructure import BandStructure

class Vsp(Prec):
    """class for reading spherical (sp) part of potential with RMT in '.vsp' file

    the potential is saved in a dict attribute, ``vsp``:

    {
        0: {(0, 0): array},
        1: {(0, 0): array},
        # ...
    }

    where array is the data, v*r.

    Latent data extraction is used here. When sp potential of atom ``ia`` is not
    requested, a tuple containing the indices of starting and end line of the data
    is placed instead of the sp data array. Upon request, the data will be extracted
    and saved.

    Args:
        pathVsp (str): the file name of vsp file
    """

    def __init__(self, pathVsp=None):
        p = pathVsp
        if p is None:
            p = get_casename() + '.vsp'
        if not os.path.isfile(p):
            raise FileNotFoundError
        
        with open(p, 'r') as h:
            self.vlines = h.readlines()
        self.vsp = {}
        self._load_data_session()
    
    def _load_data_session(self):
        """Read VLM for every atom and every LM pair

        Basically the vsp lines are first divided into sessions,
        and then each session is parsed to ``conver_lm_data`` to
        obtain the potential on each radial grid
        """
        # divide into sessions
        i = 3
        lm = None
        ia = 0
        
        while i < len(self.vlines):
            l = self.vlines[i].strip()
            if l.startswith('ATOMNUMBER'):
                if lm is not None:
                    self.vsp[ia][lm].append(i-5)
                    self.vsp[ia][lm] = tuple(self.vsp[ia][lm])
                ia = int(l.split()[-1]) - 1
                self.vsp[ia] = {}
                lm = None
            if l.startswith('VLM(R) FOR'):
                if lm is not None:
                    self.vsp[ia][lm].append(i-1)
                    self.vsp[ia][lm] = tuple(self.vsp[ia][lm])
                lm = tuple(conv_string(l, int, -3, -1))
                self.vsp[ia][lm] = [i+2,]
            i += 1
        self.vsp[ia][lm].append(i-5)
        self.vsp[ia][lm] = tuple(self.vsp[ia][lm])

    def get_vsp(self, ia):
        """Return spherical part of potential (v*r) of atom with index ``ia``

        Args:
            ia (int): index of atom, starting from 0

        Returns:
            array
        """
        lm = (0, 0)
        try:
            data = self.vsp[ia][lm]
        except KeyError:
            raise KeyError("Invalid index of atom, {} available".format(list(self.vsp.keys())))
        # load the data if the indice tuple is met
        if isinstance(data, tuple):
            st, ed = data
            self.vsp[ia][lm] = read_lm_data(self.vlines[st:ed+1])
        return self.vsp[ia][lm]


class Energy(Prec):
    """class for analysing .energy file
    """

    def __init__(self, pathEnergy=None):
        p = pathEnergy
        if p is None:
            p = get_casename() + '.energy'
        if not os.path.isfile(p):
            raise FileNotFoundError
        
        with open(p, 'r') as h:
            self.elines = h.readlines()
        if self.elines == []:
            raise IOError("Empyt energy file is parsed.")
        self.krec = 19
        self._handle_header()
        self._set_kxyz_record()
        self._read_energies()
    
    def _handle_header(self):
        """handle with the header, which consists of atomic information

        The goal is to skip the header correctly and automatically.

        It first check the scheme of float, i.e. F9.5 or F12.5
        """
        l = self.elines[0]
        # old scheme as default
        old = True
        # new scheme F12.5
        if l[3] != ".":
            if l[6] == ".":
                old = False
            else:
                raise ValueError("Invalid header format, neither F9.5 nor F12.5")
        loc = 6 - int(old) * 3
        i = 0
        while True:
            i += 1
            try:
                l = self.elines[i]
                if l[loc] != ".":
                    break
            except IndexError:
                # raise when EOF is reached but end of header is not yet reached
                raise IOError("Bad format of energy file. Check the header")
        self.nIneqAtoms = (i-1)/2
        # delete the header part
        del self.elines[:i]

    def _set_kxyz_record(self):
        """Check if extended scheme (D27.20) is used or not (ES19.12)
        """
        l = self.elines[0]
        # old scheme
        if l[2] == ".":
            self.krec = 19
        else:
            self.krec = 27
    
    def _read_energies(self):
        """Read energies in energy file
        """
        kpts = []
        weights = []
        enes = []

        i = 0
        while i < len(self.elines):
            # read the kpoint line
            l = self.elines[i]
            kpt = list(map(float, [l[self.krec * i: self.krec * (i+1)] for i in range(3)]))
            nbands, weight = conv_string(l, int, -2, -1)
            kpts.append(kpt)
            weights.append(weight)
            ene = []
            i += 1
            for j in range(nbands):
                ene.append(conv_string(self.elines[i+j], float, -1))
            enes.append(ene)
            i += nbands
        self.kpts = kpts
        self.weights = weights
        self.eigen = enes
    
    def get_nbands_min(self):
        return min([len(x) for x in self.eigen])

    def get_nbands_max(self):
        return max([len(x) for x in self.eigen])

    def load_band(self, efermi=None):
        """Return a band structure object from energy

        Args:
            efermi (float) 
        """
        raise NotImplementedError


class Scf(Prec):
    """Class for analyzing .scf file

    Args:
        pathScf (str) : the path to the case.scf file
    """

    def __init__(self, pathScf=None):
        p = pathScf
        if p is None:
            p = get_casename() + '.scf'
        if not os.path.isfile(p):
            raise FileNotFoundError
        
        with open(p, 'r') as h:
            self.scflines = h.readlines()
        # get the line indices of :ITE
        self.completed = True
        self.spin_polarized = False
        self.ites_li = find_str_matched(self.scflines, ":ITE", regex=False, full=False, ind=True)
        # check if the SCF complete
        #self._check_scf_complete()
        self.scf = {}
        self._load_iterations()
    
    def _check_scf_complete(self):
        if self.ites_li == [] or len(self.scflines) < 3:
            Verbose.print_cm_warn("SCF calculation not completed")
            self.completed = False
        if not self.scflines[-2].startswith(":ENE"):
            Verbose.print_cm_warn("SCF calculation not completed")
            self.completed = False

    def _load_iterations(self):
        for i, l in enumerate(self.scflines):
            if l.startswith(":VOL"):
                _s = self.scflines[i+2].strip()
                self.spin_polarized = _s == "SPINPOLARIZED CALCULATION"
        if self.spin_polarized:
            Verbose.print_cm_log("Spin-polarized calculation found")
            self._load_iterations_sp()
        else:
            self._load_iterations_nsp()

    def _load_iterations_nsp(self):
        infos = [("VOL", -1), ("GAP", 2), ("ENE", -1), ("FER", -1)]
        for i, iteli in enumerate(self.ites_li):
            end = len(self.scflines)
            if i < len(self.ites_li) - 1:
                end = self.ites_li[i+1]
            ite = int(self.scflines[iteli][4:7])
            oneite = self.scflines[iteli:end]
            self.scf[ite] = OrderedDict.fromkeys([x[0] for x in infos])
            try:
                for _j, l in enumerate(oneite):
                    for mark, ind in infos:
                        if l.startswith(":"+mark):
                            self.scf[ite][mark] = float(l.split()[ind])
            except ValueError:
                raise ValueError("error in loading iteration %d: %s" % (i+1, l))

    _load_iterations_sp = _load_iterations_nsp
