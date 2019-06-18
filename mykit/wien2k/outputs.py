# coding = utf-8

import os
from mykit.wien2k.utils import get_casename, read_lm_data
from mykit.core.utils import conv_string
from mykit.core.numeric import Prec
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

    def load_band(self, efermi=None):
        """Return a band structure object from energy

        Args:
            efermi (float) 
        """
        raise NotImplementedError