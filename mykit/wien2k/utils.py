# coding=utf-8

import re
import os
from shutil import which
from fnmatch import fnmatch
import numpy as np
from mykit.wien2k.constants import DEFAULT_R0, DEFAULT_R0S, DEFAULT_RMT, DEFAULT_RMTS
from mykit.core.utils import get_dirpath


def get_default_r0(elem):
    """get the default R0 for element ``elem``

    Args:
        elem (str): string representing the atomic symbol of element
            A trailing number is allowed.
    """
    e = re.sub(r"\d", '', elem)
    r0 = DEFAULT_R0S.get(e, DEFAULT_R0)
    return r0


def get_default_rmt(elem):
    """get the default RMT for element ``elem``

    Args:
        elem (str): string representing the atomic symbol of element
            A trailing number is allowed.
    """
    e = re.sub(r"\d", '', elem)
    rmt = DEFAULT_RMTS.get(e, DEFAULT_RMT)
    return rmt


def get_casename(w2kdir="."):
    """return the case name of a wien2k working directory
    
    It will first search for case.struct, if exists return the filename without extension 
    Otherwise the name of the directory will be returned
    """
    abspath = os.path.abspath(w2kdir)
    if os.path.isdir(abspath):
        for filename in os.listdir(abspath):
            if fnmatch(filename, abspath + "/*.struct"):
                case = filename.split("/")[-1][:-8]
                return case

    return os.path.basename(get_dirpath(abspath))


def find_complex_file(casename, ext):
    """Check if "{casename}.{ext}" or "{casename}.{ext}c" exists

    Args:
        casename (str) 
        ext (str) : the file extension without "c"
    """
    _path = casename + "." + ext
    if not os.path.isfile(_path):
        _path = _path + "c"
    if not os.path.isfile(_path):
        raise FileNotFoundError("Neither {casename}.in1 nor {casename}.in1c is found.")
    return _path


def get_run_lapw():
    """Return the absolute path of run_lapw executable

    Returns:
        str
    """
    return which("run_lapw")


def read_atom_info(latttype, lines):
    """Read information from struct file of a particular

    Args:
        latttype (str): the type of Bravis lattice.
            F: face-centered
            B: body-centered
        lines (list of str): the lines of file containing
            information of a particular inequivalent atom

            the first line should start with "ATOM", and
            the local rotation matrix is included
    
    Returns:
        atom, pos, r0, rmt
    """
    pos = []
    # swicth the first-atom line and mult line for convenience
    lines[0], lines[1] = lines[1], lines[0]
    mult = int(lines[0].split()[1])
    #isplit = int(lines[0].split()[3])

    for i in range(mult):
        l = lines[i + 1]
        p = list(map(float, [l[12:22], l[25:35], l[38:48]]))
        pos.append(p)
        if latttype == 'F':
            pos.append(np.add(p, [0.0, 0.5, 0.5]))
            pos.append(np.add(p, [0.5, 0.0, 0.5]))
            pos.append(np.add(p, [0.5, 0.5, 0.0]))
        elif latttype == 'B':
            pos.append(np.add(p, [0.5, 0.5, 0.5]))
    # the line including atom symbol, NPT, R0, RMT and Z
    l = lines[mult + 1]

    at = re.sub(" ", "", l[:11])
    atoms = [at,] * len(pos)
    npt = int(l[15:20])
    r0 = float(l[25:35])
    rmt = float(l[44:50])
    return atoms, pos, npt, r0, rmt


def read_symops(lines):
    """Read lines containing symmetry information in struct file

    Args:
        lines (list of str):
    
    Returns:
        dict with two keys, "rotations" and "translations"
    """
    nops = int(lines[0].split()[0])
    symops = {"rotations": np.zeros((nops, 3, 3)), "translations": np.zeros((nops, 3))}
    for i in range(nops):
        st = 4 * i + 1
        r = np.array([
            list(map(int, [lines[st + j][0:2], lines[st + j][2:4], lines[st + j][4:6]]))
            for j in range(3)
        ])
        t = np.array(list(map(float, [lines[st + j][7:] for j in range(3)])))
        symops["rotations"][i, :, :] = r
        symops["translations"][i, :] = t
    return symops
