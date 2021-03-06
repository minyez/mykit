# coding = utf-8

import os
import numpy as np

from mykit.core.log import Verbose
from mykit.core.utils import (conv_string, get_filename_wo_ext, trim_after,
                              trim_both_sides)
from mykit.wien2k.constants import (IN1_UNIT_READER_nmr, IN1_UNIT_READER_v142,
                                    IN1_UNIT_READER_v171)
from mykit.wien2k.utils import find_complex_file, get_casename

class InputError(Exception):
    pass


class In1(Verbose):
    """class for manipulating WIEN2k in1 file

    TODO:
        Read data in and after KPOINT UNIT line
    
    Args:
        casename (str)
        switch (str)
        efermi (float) : Fermi energy
        rkmax (float) : RKmax
        lmax (int) : the maximum angular quantum of APW basis in MT region
        lnsmax (int)
        cmplx (bool)
        elparams (dict): the various information of exceptions of each atom
    """

    def __init__(self, casename, switch, efermi, rkmax, lmax, lnsmax, cmplx,
                 unit, emin, emax, nbands, *elparams, **kwargs):
        self.casename = casename
        self.switch = switch
        self.efermi = efermi
        self.rkmax = rkmax
        self.lmax = lmax
        self.lnsmax = lnsmax
        self.cmplx = cmplx
        self.unit = unit
        self.emin = emin
        self.emax = emax
        self.nbands = nbands
        self.elparams = elparams

    def add_exception(self, atomId, l, e, search=0.000, cont=True, apw=1):
        """Add one exception for a particular atom and l channel

        Args:
            atomId (int): index of inquivalent atom
            l (int): the angular momentum
            e (float): the linearization energy to set
            search (float): step size for energy search
            cont (bool)
            apw (0, 1)
        """
        # add to existing item only
        try:
            elparam = self.elparams[atomId]
        except IndexError:
            raise IndexError("atom index %d not available" % atomId)
        assert isinstance(cont, bool)
        assert apw in [0, 1]
        contStr = {True: "CONT", False: "STOP"}[cont]
        elparam["Ndiff"] += 1
        if not l in elparam["exceptions"]:
            elparam["exceptions"][l] = []
        elparam["exceptions"][l].append([e, search, contStr, apw])
    
    def remove_last_exception(self, atomId, l):
        """Remove the last linearization energy of atomId in channel l

        Args:
            atomId (int): index of inquivalent atom
            l (int): the angular momentum
        """
        try:
            elparam = self.elparams[atomId]
        except IndexError:
            raise IndexError("atom index %d not available" % atomId)
        if l in elparam["exceptions"]:
            elparam["Ndiff"] -= 1
            del elparam["exceptions"][l][-1]
    
    @property
    def natoms(self):
        """int. number of inequivalent atoms
        """
        return len(self.elparams)

    def get_exceptions(self, atomId, l=None, eonly=False):
        """Get exceptions of a particular atom

        Args:
            atomId (int)
            l (int): if specified, the exceptions in particular l channel is returned
            eonly (bool): get the linearization energy only
        
        Returns:
            dict if l is None. Otherwise list.
        """
        assert isinstance(atomId, int)
        try:
            elparam = self.elparams[atomId]["exceptions"]
            if eonly:
                es = {}
                for l in elparam.keys():
                    e = [x[0] for x in elparam[l]]
                    es[l] = e
                elparam = es
        except IndexError:
            return None
        if l is not None:
            elparam = elparam.get(l, None)
        return elparam
    
    def __str__(self):
        s = ["%5s  EF=%12.11f" % (self.switch, self.efermi), ]
        s += ["%4.1f %8d %4d" % (self.rkmax, self.lmax, self.lnsmax)]
        for _ia, elparams in enumerate(self.elparams):
            s.append(_write_el_block(elparams))
        last = "K-VECTORS FROM UNIT:%1d %6.1f" % (self.unit, self.emin)
        if self.emax > 100:
            last += "%10.1E" % self.emax
        else:
            last += "%10.1f" % self.emax
        if self.nbands is not None:
            last += "%6d" % self.nbands
        s.append(last)
        return '\n'.join(s)

    def write(self, pathIn=None, backup=False, suffix="_bak"):
        """Write to in1 file
        """
        path = pathIn
        if path is None:
            path = self.casename + 'in1' + 'c' * int(self.cmplx)
        if os.path.isdir(path):
            raise IOError("Trying to write in1 as an existing directory")
        if os.path.isfile(path) and backup:
            bakpath = path + suffix.strip()
            os.rename(path, bakpath)
        with open(path, 'w') as f:
            print(self.__str__(), file=f)

    @classmethod
    def read_from_file(cls, pathIn1=None):
        """Return In1 instance by reading an exisiting 

        Args:
            filePath (str): the file name
        """
        cmplx = False
        if pathIn1 is None:
            casename = get_casename()
            path, cmplx = find_complex_file(casename, "in1")
        else:
            path = pathIn1
            if path.endswith("c"):
                cmplx = True
            casename = get_filename_wo_ext(path)

        with open(path, "r") as h:
            w2klines = h.readlines()
        switch = w2klines[0][:5]
        efermi = float(trim_both_sides(w2klines[0], r"=", r"\("))
        _params = trim_after(w2klines[1], r"\(").split()
        rkmax = float(_params[0])
        lmax, lnsmax = tuple(map(int, _params[1:3]))

        elparams = []
        i = 2
        _flagEnergy = True
        _flagInExcept = False
        while i < len(w2klines):
            line = trim_after(w2klines[i], r"\(")
            if line.startswith("K-VECTORS FROM UNIT"):
                _flagEnergy = False
                try:
                    # wien2k 17.1
                    data = IN1_UNIT_READER_v171.read(line)
                    unit, emin, emax, nbands = data
                    de = emax - emin
                except ValueError:
                    try:
                        # wien2k 14.2
                        data = IN1_UNIT_READER_v142.read(line)
                        unit, emin, de, nbands = data
                        emax = efermi + de
                    except ValueError:
                        # meet NMR in1
                        data = IN1_UNIT_READER_nmr.read(line)
                        unit, emin, emax = data
                        nbands = None
            else:
                words = line.split()
                if _flagEnergy and len(words) == 3:
                    ndiff = int(words[1])
                    atomEl = _read_el_block(w2klines[i : i + ndiff + 1])
                    elparams.append(atomEl)
                    i += ndiff
            i += 1
        return cls(casename, switch, efermi, rkmax, lmax, lnsmax, cmplx,
                   unit, emin, emax, nbands, *elparams)


def _read_el_block(elBlock):
    """
    Args:
        elBlock (list or tuple):  strings containing el information

    Returns
        dict
    """
    try:
        assert isinstance(elBlock, (list, tuple))
        for i in elBlock:
            assert isinstance(i, str)
    except AssertionError:
        raise InputError("elBlock should be a list of strings")
    elParams = {}
    globParams = elBlock[0].split()
    elParams["Etrial"] = float(globParams[0])
    ndiff = int(globParams[1])
    if len(elBlock) != ndiff + 1:
        raise InputError(
            "Inconsistent El block: need {}, parsed {}".format(ndiff, len(elBlock) - 1)
        )
    elParams["Ndiff"] = ndiff
    elParams["Napw"] = int(globParams[2])
    exceptions = {}
    for j in range(ndiff):
        # l, e, eIncre, cont, apw = EL_READER.read(elBlock[1 + j])
        s = elBlock[1+j]
        l, apw = conv_string(s, int, 0, -1)
        e, eIncre = conv_string(s, float, 1, 2)
        cont = s.split()[-2]
        if not l in exceptions:
            exceptions[l] = []
        exceptions[l].append([e, eIncre, cont, apw])
    elParams["exceptions"] = exceptions
    return elParams


def _write_el_block(elParams):
    """
    Args:
        elParams (dict)

    Returns
        str
    """
    ep = elParams
    ret = ["%6.2f %4d %2d" % (ep["Etrial"], ep["Ndiff"], ep["Napw"])]
    # sort by l
    exceptions = sorted(list(ep["exceptions"].items()), key=lambda x: x[0])
    for l, els in exceptions:
        for el in els:
            ret.append("%2d %9.5f %9.5f %4s %1d" % (l, *el))
    return "\n".join(ret)

