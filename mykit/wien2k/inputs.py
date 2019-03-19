# coding = utf-8

import os

from fortranformat._exceptions import InvalidFormat

from mykit.core.utils import (get_cwd_name, get_filename_wo_ext, trim_after,
                              trim_both_sides)
from mykit.wien2k.constants import EXCEPTION_READER, EXCEPTION_WRITER


class InputError(Exception):
    pass


class In1:
    '''class for manipulating WIEN2k in1 file

    TODO:
        Read data in and after KPOINT UNIT line
    '''
    def __init__(self, casename, switch, efermi, \
                       rkmax, lmax, lnsmax, \
                    #    unit, emin, emax, nband, \
                       *elparams):
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
        pass

    def write(self):
        pass

    @classmethod
    def read_from_file(cls, filePath):
        '''Return In1 instance by reading an exisiting 

        Args:
            filePath (str): the file name
        '''
        if filePath is None:
            casename = get_cwd_name()
            _path = casename + '.in1'
            if not os.path.isfile(_path):
                _path = casename + '.in1c'
            if not os.path.isfile(_path):
                raise InputError("Neither {casename}.in1 nor {casename}.in1c is found.")
        else:
            _path = filePath
            casename = get_filename_wo_ext(_path)
        
        with open(_path, 'r') as h:
            w2klines = h.readlines()
        switch = w2klines[0][:5]
        efermi = float(trim_both_sides(w2klines[0], r"=", r"\("))
        params = trim_after(w2klines[1], r"\(").split()
        rkmax = float(params[0])
        lmax, lnsmax = tuple(map(int, params[1:3]))

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
                    atomExcept = {}
                    atomExcept["Etrial"] = float(words[0])
                    ndiff = int(words[1])
                    atomExcept["Ndiff"] = ndiff
                    atomExcept["Napw"] = int(words[2])
                    exceptions = {}
                    for j in range(ndiff):
                        l, e, eIncre, cont, apw = EXCEPTION_READER.read(w2klines[i+1+j])
                        if not l in exceptions:
                            exceptions[l] = []
                        exceptions[l].append((e, eIncre, cont, apw))
                    atomExcept["exceptions"] = exceptions
                    elparams.append(atomExcept)
                    i += ndiff
            i += 1
        return cls(casename, switch, efermi, rkmax, lmax, lnsmax, \
                    # unit, emin, emax, nband, \
                    *elparams)


def read_el_exception_block(elBlock):
    pass


def write_el_exception_block(elParams):
    pass
