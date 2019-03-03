# -*- coding: utf-8 -*-
'''include class for KPOINTS
'''

import os
from mykit.core.kmesh import kmesh_control
from mykit.core.log import verbose

class KpointsError(Exception):
    pass

# ! kmesh_control instance is used as an attribute,
# ! since vasp does not have any tag for it
class kpoints(verbose):
    '''class for manipulation KPOINTS file.

    Note:
        Only automatic generation is supported.
        And it does not support modification of arguments

    TODO:
        Line mode

    Args:
        comment (str): the comment for KPOINTS
        kgrids (3-member list): the number of grid on each 
        kdense (float): the density of k grid point. 
            If set to non-zero value and kmode is "G", "M" or "A", the kmode will automatically be set as "A".
            If kmode is set as "L" and kpath specified, it is the number of kpoints along two special points.
        kshift (3-member list): shift off center
        kmode ("G", "M", "A", "L"): the mode of KPOINTS
            "G": Gamma-centered
            "M": Monkhorse-Pack
            "A": fully automatic, i.e. the grid will be decided by kdense
            "L": line mode for bandstructure. kpath should be parse as well. Reciprocal vector unit is always used.
        kpath (nx3 array): the special points along the line.
        kpoints (nx4 array): the explicit kpoints to parse. the fourth column is the kpoint weight.
    '''

    def __init__(self, comment=None, kgrids=None, kdense=None, kshift=None, kmode="G", kpath=None, kpoints=None):
        self.comment = ''
        # raise when no usable parameter is input
        if not (kgrids or kdense or kpath or kpoints):
            raise KpointsError("enter at least one of kgrids, kdense, kpath and kpoints")
        self._control = kmesh_control("mykit")
        self._relatePoscar = None
        # Search for POSCAR when kdense is specified
        if kdense != None:
            from mykit.vasp.poscar import poscar, PoscarError
            try:
                self._poscar = poscar.read_from_file("POSCAR")
            except PoscarError:
                raise KpointsError("POSCAR not found while kmesh density is set.")

    def __print(self, fp):
        '''

        Args:
            fp (file)
        '''
        pass

    def print(self):
        '''Preview the KPOINTS output
        '''
        from sys import stdout
        self.__print(stdout)
        
    def write(self, pathKpoints="KPOINTS", backup=False, suffix="_bak"):
        '''Write KPOINTS to path
        '''
        _name = pathKpoints
        try: 
            assert not os.path.isdir(_name)
        except AssertionError:
            raise KpointsError("The path to write KPOINTS is a directory.")
        if os.path.isfile(_name) and backup:
            _bakname = _name + suffix.strip()
            os.rename(_name, _bakname)
        with open(_name, 'w') as f:
            self.__print(f)

    @classmethod
    def read_from_file(cls, pathKpoints="KPOINTS"):
        pass