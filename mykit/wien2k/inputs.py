# coding = utf-8

import os

from fortranformat._exceptions import InvalidFormat

from mykit.core.log import Verbose
from mykit.core.utils import get_filename_wo_ext, trim_after, trim_both_sides
from mykit.wien2k.utils import get_casename, find_complex_file
from mykit.wien2k.constants import EL_READER, EL_WRITER


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
        elparams (dict): the various information of exceptions of each atom
    """

    def __init__(
        self,
        casename,
        switch,
        efermi,
        rkmax,
        lmax,
        lnsmax,
        #    unit, emin, emax, nband, \
        *elparams,
        **kwargs
    ):
        self.casename = casename
        self.switch = switch
        self.efermi = efermi
        self.rkmax = rkmax
        self.lmax = lmax
        self.lnsmax = lnsmax
        # self.unit = unit
        # self.emin = emin
        # self.emax = emax
        # self.nband = nband
        self.elparams = elparams

    def add_exception(self, atomId, l, e, search=0.000, cont=True, lapw=0):
        assert isinstance(cont, bool)
        assert isinstance(atomId, int)
        assert lapw in [0, 1]
        contStr = {True: "CONT", False: "STOP"}[cont]
        # add to existing item only
        elparam = self.elparams[atomId]
        elparam["Ndiff"] += 1
        if not l in elparam["exceptions"]:
            elparam["exceptions"][l] = []
        elparam["exceptions"][l].append([e, search, contStr, lapw])

    def write(self):
        # TODO
        raise NotImplementedError

    @classmethod
    def read_from_file(cls, filePath=None):
        """Return In1 instance by reading an exisiting 

        Args:
            filePath (str): the file name
        """
        if filePath is None:
            casename = get_casename()
            _path = find_complex_file(casename, "in1")
        else:
            _path = filePath
            casename = get_filename_wo_ext(_path)

        with open(_path, "r") as h:
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
                # try:
                #     # wien2k 17.1
                #     data = IN1_UNIT_READER_v171.read(line)
                #     unit, emin, emax, nband = data
                #     de = emax - emin
                # except InvalidFormat:
                #     # wien2k 14.2
                #     data = IN1_UNIT_READER_v142.read(line)
                #     unit, emin, de, nband = data
                #     emax = efermi + de
            else:
                words = line.split()
                if _flagEnergy and len(words) == 3:
                    ndiff = int(words[1])
                    atomEl = read_el_block(w2klines[i : i + ndiff + 1])
                    elparams.append(atomEl)
                    i += ndiff
            i += 1
        return cls(
            casename,
            switch,
            efermi,
            rkmax,
            lmax,
            lnsmax,
            # unit, emin, emax, nband, \
            *elparams
        )


def read_el_block(elBlock):
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
        l, e, eIncre, cont, apw = EL_READER.read(elBlock[1 + j])
        if not l in exceptions:
            exceptions[l] = []
        exceptions[l].append([e, eIncre, cont, apw])
    elParams["exceptions"] = exceptions
    return elParams


def write_el_block(elParams):
    """
    Args:
        elParams (dict)

    Returns
        str
    """
    _ret = []
    # TODO
    return "\n".join(_ret)
