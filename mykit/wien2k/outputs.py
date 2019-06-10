# coding = utf-8

import os
from mykit.wien2k.utils import get_casename, read_lm_data
from mykit.core.utils import conv_string
from mykit.core.numeric import Prec

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
