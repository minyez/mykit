# coding=utf-8

import re
import os
from shutil import which
from fnmatch import fnmatch
import numpy as np
from mykit.wien2k.constants import DEFAULT_R0, DEFAULT_R0S, DEFAULT_RMT, DEFAULT_RMTS, CFEIN2
from mykit.core.utils import get_dirpath
from mykit.core.elements import NUCLEAR_CHARGE
from mykit.core.numeric import Prec


def get_default_r0(elem):
    """get the default R0 for element ``elem``

    Args:
        elem (str): string representing the atomic symbol of element
            A trailing number is allowed.
    """
    e = re.sub(r"\d", '', elem)
    return DEFAULT_R0S.get(e, DEFAULT_R0)


def get_default_rmt(elem):
    """get the default RMT for element ``elem``

    Args:
        elem (str): string representing the atomic symbol of element
            A trailing number is allowed.
    """
    e = re.sub(r"\d", '', elem)
    return DEFAULT_RMTS.get(e, DEFAULT_RMT)


def get_z(elem):
    """Get nuclear charge Z from element ``elem``

    Args:
        elem (str): string representing the atomic symbol of element
            A trailing number is allowed.
    """
    e = re.sub(r"\d", '', elem)
    z = NUCLEAR_CHARGE.get(e, None)
    if z is None:
        raise ValueError("Invalid element identifier: %s" % elem)
    return float(z)


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
    path = casename + "." + ext
    c = False
    if not os.path.isfile(path):
        path += "c"
        c = True
    if not os.path.isfile(path):
        raise FileNotFoundError("Neither {casename}.{ext} nor {casename}.{ext}c is found.")
    return path, c


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


def read_lm_data(lines, lenfield=19):
    """Read the potential data of a particular LM channel in vsp and vns file

    Note that the header containing 'L=,M=' should be excluded. Only numbers and
    empty lines are allowed. Unit conversion from Rydberg to Hartree is done.

    Args:
        lines (list of str): file lines containing data of potential
        lenfield (int): the length of field that one number occupies, 19 as default

    Returns:
        list, potential in Hartree unit
    """
    l = ''.join([x.strip() for x in lines])
    ll = len(l)
    if ll%lenfield != 0:
        raise IOError("Bad input, total length %d not divided by %d" % (ll, lenfield))
    n = int(ll/lenfield)
    data = [float(l[lenfield*i:lenfield*(i+1)])/2.0 for i in range(n)]
    return data


def solve_rde(r, vr, E, l, Z, normalize=False):
    """Python adaption of OUTWINB subroutine for outward integration
    of radial Dirac equation. Scalar relativistic effect is always
    included.

    Args
        r (array-like): radial grid points
        vr (array-like): potential times r, V*r, in a.u.
        E (float): energy in Hatree
        l (int): angular momentum
        Z (float): nuclear charge

    Returns
        array, array, float, float, int: 
            u*r(large), u*r(small), u and du/dr at boundary, and number of nodes
    """
    ngrids = len(r)
    nodes = 0
    ERy = E * 2.0
    assert len(vr) == ngrids
    assert Z > 0.9
    assert isinstance(E, float)
    assert isinstance(l, (float, int))

    ur_large = np.zeros(ngrids, dtype=Prec._dtype)
    ur_small = np.zeros(ngrids, dtype=Prec._dtype)
    r = np.array(r, dtype=Prec._dtype)
    v = 2.0 * np.array(vr, dtype=Prec._dtype)/r
    dx = np.log(r[-1]/r[0]) / (ngrids-1)
    c = 2.0e0 * 137.0359895
    z = float(Z) * 2
    zdc = z/c

    llp1 = float(l*(l+1))
    s = np.sqrt(llp1 + 1.0 - zdc**2)

    dgf = np.zeros((2, 3))
    ur_large[0:3] = np.power(r[0:3], s) * 1.0
    ur_small[0:3] = np.power(r[0:3], s) * (s-1) / zdc
    dgf[0, :] = dx * ur_large[0:3] * s
    dgf[1, :] = dx * ur_small[0:3] * s
    dg1 = dgf[0, 0]
    dg2 = dgf[0, 1]
    dg3 = dgf[0, 2]
    df1 = dgf[1, 0]
    df2 = dgf[1, 1]
    df3 = dgf[1, 2]

    for i in range(3, ngrids):
        rdx = r[i] * dx
        phi = (ERy - v[i]) * rdx / c
        U = rdx * c + phi
        Y = - llp1 * dx * dx / U + phi
        det = 64.0/9.0 - dx * dx + U * Y
        B11 = ur_large[i-1] * 8.0/3.0 + dg1 / 9.0 - 5.0 * dg2 / 9.0 + 19.0 * dg3 / 9.0
        B22 = ur_small[i-1] * 8.0/3.0 + df1 / 9.0 - 5.0 * df2 / 9.0 + 19.0 * df3 / 9.0
        ur_large[i] = (B11*(8.0/3.0+dx) + B22 * U)/det
        ur_small[i] = (B22*(8.0/3.0-dx) - B11 * Y)/det
        if ur_large[i-1]*ur_large[i] < 0:
            nodes += 1
        dg1, dg2, dg3 = dg2, dg3, U*ur_small[i] + dx * ur_large[i]
        df1, df2, df3 = df2, df3, - Y*ur_large[i] - dx * ur_small[i]

    ur_small *= c/2.0

    u = ur_large[-1]/r[-1]
    dudr = (dg3/dx/r[-1] - u)/r[-1]
    if normalize:
        norm = np.sqrt(logr_int(r, ur_large, ur_small, ur_large, ur_small))
        ur_large /= norm
        ur_small /= norm
        u /= norm
        dudr /= norm
    return ur_large, ur_small, u, dudr, nodes


def solve_ene_deriv(r, vr, E, l, Z):
    """Compute the energy derivative of wavefunction at E
    using finite difference
    """
    dE = 1.0e-3
    urlp, ursp, up, dup, _ = solve_rde(r, vr, E+dE, l, Z, normalize=True)
    urln, ursn, un, dun, _ = solve_rde(r, vr, E-dE, l, Z, normalize=True)
    n = len(r)
    urlde = (urlp - urln) / 2.0 / dE
    ursde = (ursp - ursn) / 2.0 / dE
    ude = (up - un) / 2.0 / dE
    dude = (dup - dun) / 2.0 / dE
    # normalize to u at E
    url, urs, u, du, _ = solve_rde(r, vr, E, l, Z, normalize=True)
    ovlp = logr_int(r, urlde, ursde, url, urs)
    urlde -= ovlp * url
    ursde -= ovlp * urs
    ude -= ovlp * u
    dude -= ovlp * du
    nodes = 0
    for i in range(1, n):
        if urlde[i-1] * urlde[i] < 0.0:
            nodes += 1
        if urlde[i-1] == 0.0:
            nodes -= 1
    return urlde, ursde, ude, dude, nodes


def logr_int(r, a, b, x, y):
    """Calculate the sum of inner product, (a, x) + (b, y)
    on a logarithmic radial grid point.

    This is adapted from rint13 in WIEN2k
    """
    def _func(n):
        return r[n] * (a[n] * x[n] + CFEIN2 * b[n] * y[n])

    ngrids = len(r)
    assert ngrids == len(a) == len(b) == len(x) == len(y)
    dx = np.log(r[1]/r[0])
    j = 3 - ngrids % 2
    j1 = j - 1
    z4 = 0.0
    z2 = 0.0
    while True:
        z4 += _func(j-1)
        j += 1
        if j >= ngrids:
            break
        z2 += _func(j-1)
        j += 1
    P1 = _func(0)
    P2 = _func(j1-1)
    inte = 2.0 * z2 + 4.0 * z4 + P2 + _func(-1)
    inte = (dx*inte + P1)/3.0
    if j1 > 1:
        inte += 0.5 * dx * (P1+P2)
    return inte
