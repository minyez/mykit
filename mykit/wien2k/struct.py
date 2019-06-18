# coding = utf-8
"""Class that manipulates WIEN2k main input struct file 
"""
import numpy as np
from mykit.core.cell import Cell
#from mykit.core.symmetry import get_spacegroup, get_sym_operations, standardize
from mykit.core.utils import (
    get_latt_vecs_from_latt_consts,
    get_all_atoms_from_sym_ops,
)
from mykit.wien2k.constants import STRUCT_LATT_PARAM_READER, DEFAULT_NPT
from mykit.wien2k.utils import read_atom_info, read_symops, get_default_r0, get_default_rmt, get_z, get_casename


class StructError(Exception):
    pass


class Struct(Cell):
    """Class to manipulate WIEN2k master input case.struct file

    Args:
        latt, atoms, pos, kwargs: see docstring of Cell
        r0 (dict)
        rmt (dict)
    """

    def __init__(self, latt, atoms, pos, npt=None, r0=None, rmt=None, **kwargs):
        super(Struct, self).__init__(latt, atoms, pos, **kwargs)
        self.r0 = r0
        self.rmt = rmt
        self.npt = npt
        self._set_r0()
        self._set_rmt()
        self._set_npt()

    def _set_r0(self):
        if self.r0 is None:
            r0 = {}
            for a in self.atomTypes:
                r0[a] = get_default_r0(a)
            self.r0 = r0
        else:
            assert isinstance(self.r0, dict)
            for a in self.atomTypes:
                if a not in self.r0:
                    self.r0[a] = get_default_r0(a)

    def _set_rmt(self):
        if self.rmt is None:
            rmt = {}
            for a in self.atomTypes:
                rmt[a] = get_default_rmt(a)
            self.rmt = rmt
        else:
            assert isinstance(self.rmt, dict)
            for a in self.atomTypes:
                if a not in self.rmt:
                    self.rmt[a] = get_default_rmt(a)

    def _set_npt(self):
        if self.npt is None:
            npt = {}
            for a in self.atomTypes:
                npt[a] = DEFAULT_NPT
            self.npt = npt
        else:
            assert isinstance(self.npt, dict)
            for a in self.atomTypes:
                if a not in self.npt:
                    self.npt[a] = DEFAULT_NPT

    def __str__(self):
        raise NotImplementedError

    def _get_atom_from_iden(self, iden):
        if isinstance(iden, int):
            at = self.atomTypes[iden]
        elif isinstance(iden, str):
            if iden in self.atomTypes:
                raise ValueError("atom %s is not found" % iden)
            at = iden
        else:
            raise TypeError("iden should be int or str")
        return at

    def get_radial(self, iden):
        """return the radial grid points of particular muffin-tin

        Args:
            iden (int or str): the identifier for inequivalent atom
                int for the index and str for the name.
        """
        at = self._get_atom_from_iden(iden)
        r0 = self.r0[at]
        rmt = self.rmt[at]
        npt = self.npt[at]
        return r0 * np.logspace(0, np.log(rmt/r0), num=npt, base=np.e)
    
    def get_z(self, iden):
        """return the nuclear charge of atom

        Args:
            iden (int or str): the identifier for inequivalent atom
                int for the index and str for the name.
        """
        at = self._get_atom_from_iden(iden)
        return get_z(at)
    
    @classmethod
    def read_from_file(cls, pathStruct=None):
        """Read from an existing struct file

        Args:
            pathStruct (str): path to the file to read as WIEN2k struct
        """
        p = pathStruct
        if p is None:
            p = get_casename() + ".struct"

        with open(p, "r") as h:
            slines = h.readlines()

        kwargs = {"coordSys": "D", "unit": "au"}
        kwargs["comment"] = slines[0].strip()
        nIneqAtoms = int(slines[1].split()[2])
        latttype = slines[1].split()[0]
        # TODO: mode of calculation
        latt = get_latt_vecs_from_latt_consts(*STRUCT_LATT_PARAM_READER.read(slines[3]))

        ineqAtoms = []
        ineqPos = []
        npt = {}
        r0 = {}
        rmt = {}

        atomBlocks = []
        newAtom = True
        # divide lines into block of atoms and symmetry operations
        for i, line in enumerate(slines):
            if i < 4:
                continue
            l = line.strip()
            if l.startswith("ATOM") and newAtom:
                s = i
                newAtom = False
            if l.startswith("LOCAL ROT MATRIX"):
                atomBlocks.append(tuple([s, i + 2]))
                newAtom = True
            if l.endswith("SYMMETRY OPERATIONS"):
                symopsSt = i
                break
        assert len(atomBlocks) == nIneqAtoms

        for s, l in atomBlocks:
            at, pos, anpt, ar0, armt = read_atom_info(latttype, slines[s : l + 1])
            ineqAtoms.extend(at)
            ineqPos.extend(pos)
            npt[at[0]] = anpt
            r0[at[0]] = ar0
            rmt[at[0]] = armt
        symops = read_symops(slines[symopsSt:])
        atoms, pos = get_all_atoms_from_sym_ops(ineqAtoms, ineqPos, symops)
        return cls(latt, atoms, pos, npt, r0, rmt, **kwargs)


