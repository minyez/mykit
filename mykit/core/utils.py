# -*- coding: utf-8 -*-
'''this module defines some common used utilities'''


import os
import re
import subprocess as sp
from collections import OrderedDict
from collections.abc import Iterable
from shutil import rmtree
from sys import stdout

import CifFile
import numpy as np


def get_dirpath(filePath):
    '''get the name of directory with filePath

    Args:
        filePath (str): the string of the path of file

    Returns:
        str: the absolute path of parent directory, if filepath represents a file,
            otherwise the path of the directory
    '''

    _path = os.path.abspath(filePath)
    if os.path.isdir(_path):
        _path = _path + '/'
    return os.path.dirname(_path)


def get_file_ext(filePath):
    '''Return the extension name of filePath

    If filePath is a existing directory, None will be returned
    If the path have no characters after "." or have no ".", 
    an empty string will be returned.

    Args:
        filePath (str): the path of the file
    '''
    if os.path.isdir(filePath):
        return None
    base = os.path.basename(os.path.abspath(filePath))
    return os.path.splitext(base)[1][1:]


def get_filename_wo_ext(filePath):
    '''Get the filename without extension

    Args:
        filePath (str): the path of file
    '''
    fnExt = os.path.basename(os.path.abspath(filePath))
    return os.path.splitext(fnExt)[0]


def get_cwd_name():
    '''Get the name of current working directory
    '''
    return os.path.basename(os.getcwd())


def get_matched_files(dirPath='.', regex=None):
    '''Get the abspath of the files whose name matches a regex

    Only files will be returned, and directories are excluded.

    Args:
        dirPath (str): the directory to search
        regex (regex): the regular expression to match the filename

    Returns:
        tuple of strings
    '''
    # check the exisitence of path
    fns = []
    _absDir = os.path.abspath(dirPath)
    if os.path.isdir(_absDir):
        for i in os.listdir(_absDir):
            if regex != None:
                if not re.match(regex, i):
                    continue
            _fpath = os.path.join(_absDir, i)
            if os.path.isfile(_fpath):
                fns.append(_fpath)
    return tuple(fns)


# def common_io_checkdir(dirname=None, create=True):
#     '''
#     check if dirname exists, or create it
#     return: the full path of target directory
#     '''
#     dirname = dirname.strip()

#     if (dirname is None or dirname.strip() == ""):
#         dirname = os.getcwd()
#     elif (not os.path.exists(dirname)) and create:
#         os.mkdir(dirname.strip())
#     return dirname


# def common_io_cleandir(dirname=None):
#     '''
#     check if dirname exists and is empty, or create it
#     and make it contain no files
#     return: the full path of target directory
#     '''
#     if (dirname is None or dirname.strip() == ""):
#         dirname = os.getcwd()
#     elif (not os.path.exists(dirname)):
#         os.mkdir(dirname)
#     elif (os.path.exists(dirname)):
#         rmtree(dirname)
#         os.mkdir(dirname)
#     return dirname


def trim_after(string, regex, include_pattern=False):
    '''Trim a string after the first match of regex.

    If fail to match any pattern, the original string is returned

    The matched pattern is trimed as well.

    Args:
        string (str): the string to trim
        regex (regex): the regex to match
        include_pattern (bool): if the matched pattern is included
        in the return string
    '''
    _search = re.search(regex, string)
    if _search != None:
        if include_pattern:
            return string[:_search.end()]
        else:
            return string[:_search.start()]
    return string


def trim_before(string, regex, include_pattern=False):
    '''Trim a string from the beginning to the first match of regex.

    If fail to match any pattern, the original string is returned.

    Args:
        string (str): the string to trim
        regex (regex): the regex to match
        include_pattern (bool): if the matched pattern is included
        in the return string
    '''
    _search = re.search(regex, string)
    if _search != None:
        if include_pattern:
            return string[_search.start():]
        else:
            return string[_search.end():]
    return string


def trim_both_sides(string, regex_left, regex_right, include_pattern=False):
    '''Trim a string from both sides.

    Basically it first tries to match regex_left, trim the characters on the left
    of the matched pattern, then match regex_right and trim the characters after.

    Args:
        regex_left (regex):
        regex_right (regex):
        include_pattern (bool): if the matched pattern is included
        in the return string
    '''
    _trimed = trim_before(string, regex_left, include_pattern=include_pattern)
    _trimed = trim_after(_trimed, regex_right, include_pattern=include_pattern)
    return _trimed


def check_duplicates_in_tag_tuple(tagtuple):
    '''Check if there is duplicate in a tag tuple, case sensitive

    Args:
        tagTuple (tuple) : the tag tuple to check
    '''
    _dup = -1
    for _i, _k in enumerate(tagtuple):
        if _k in tagtuple[:_i]:
            _dup = _i
            break
    return _dup


# def data_normalization(data, scale=1.0, normByPeak=True):
#     '''Normalize the 1D data.

#     Args:
#         data (iterable): the container of 1D data
#         normByPeak (bool) : when set True, the normalization factor will be
#             the peak absolute value. Otherwise, the sum of absolute values
#             will be used as normalization factor.

#     TODO:
#         Generalize the normalization

#     Returns:
#         numpy array, the normalized data
#     '''
#     import numpy as np
#     assert len(np.shape(data)) == 1
#     assert isinstance(normByPeak, bool)

#     _a = []
#     try:
#         _a = np.array(data, dtype="float64")
#     except ValueError:
#         raise ValueError("the data cannot be converted.")

#     _sum = np.sum(np.abs(_a)) / scale
#     _max = np.max(np.abs(_a)) / scale
#     if normByPeak:
#         return _a / _max
#     return _a / _sum


# def find_data_extreme(data):
#     '''Find the point at which the data reaches extrema

#     TODO:
#         Generalize to 2D and 3D coordinate

#     Returns:
#         dict, with two keys, "min" and "max".
#         Either key has a 2-member tuple with its first the min/max value
#         and second the coordinate where it reaches the extreme
#     '''
#     pass


def find_vol_dirs(path='.', vdPat=None):
    '''Find names of directories corresponding to calculation with lattice of different volumes

    Args:
        path (str): the path to search directories within. Default is CWD.
        vdPat (regex): the pattern of the names of volume directories
            If not specified, use "V_x.xx" where x is 0-9

    Returns:
        list of strings
    '''
    pat = vdPat
    if pat is None:
        pat = r'^V_\d.\d+'
    _dirs = []
    for _d in os.listdir(path):
        if re.match(pat, _d):
            _dirs.append(_d)

    def __sort_vol(dirstr):
        return float(dirstr.split('_')[1])
    if vdPat is None:
        _dirs = sorted(_dirs, key=__sort_vol)
    return _dirs


def conv_string(string, conv2, *indices, sep=None, strips=None):
    '''
    Split the string and convert substrings to a specified type.

    Args:
        string (str): the string from which to convert value
        conv2: the type to which the substring will be converted
            support ``str``, ``int``, ``float``, ``bool``
        indices (int): if specified, the substring with indices in the splitted string lists will be converted.
            otherwise, all substring will be converted.
        sep (regex): the separators used to split the string.
        strips (str): extra strings to strip for each substring before conversion

    Returns:
        ``conv2``, or list of ``conv2`` type
    '''
    assert conv2 in [str, int, float, bool]
    str_tmp = string.strip()
    if sep is not None:
        str_list = re.split(sep, str_tmp)
    else:
        str_list = str_tmp.split()
    if strips is None:
        str_list = [x.strip() for x in str_list]
    else:
        str_list = [x.strip(' '+strips) for x in str_list]

    # need to convert to float first for converting to integer
    if conv2 is int:
        def convfunc(x): return int(float(x))
    elif conv2 is bool:
        def convfunc(x): return {'TRUE': True, 'T': True, '.TRUE.': True, '.T.': True,
                                 'FALSE': True, 'F': True, '.FALSE.': True, '.F.': False, }.get(x.upper(), None)
    else:
        convfunc = conv2

    if len(indices) == 0:
        return list(map(convfunc, str_list))
    elif len(indices) == 1:
        return convfunc(str_list[indices[0]])
    else:
        conv_strs = [str_list[i] for i in indices]
        return list(map(convfunc, conv_strs))


# def common_ss_conv(string, i, conv2, sep=None):
#     '''
#     Split the string and convert a single substring to a specified type.

#     Args:
#         string (str): the string from which to convert value
#         i (int): the substring index in the list to be converted after splitting by sep
#         conv2: the type to which the substring will be converted
#         sep (regex): the separators used to split the string.
#     '''

#     str_tmp = string.strip()
#     if sep is not None:
#         str_list = re.split(r'[%s]'%sep, str_tmp)
#     #    print(str_list)
#     else:
#         str_list = str_tmp.split()

#     # try:
#         # return conv2(str_list[i])
#     # except ValueError:
#     #     return conv2(float(str_list[i]))
#     return conv2(str_list[i])


def get_first_last_line(filePath, encoding=stdout.encoding):
    '''Return the first and the last lines of file

    The existence of filePath should be check beforehand.

    Args:
        filePath (str): the path of the file
        encoding (str): the encoding of the file. Default stdout.encoding

    Returns
        two strings (unstripped)
    '''
    with open(filePath, "rb") as f:
        first = f.readline()        # Read the first line.
        f.seek(-2, os.SEEK_END)     # Jump to the second last byte.
        while f.read(1) != b"\n":   # Until EOL is found...
            # ...jump back the read byte plus one more.
            f.seek(-2, os.SEEK_CUR)
        last = f.readline()         # Read last line.
    # encode string
    return str(first, encoding), str(last, encoding)


def get_str_indices(container, string):
    '''Return the indices of ``string`` in a list or tuple``container``

    Args:
        container (list or tuple): container of strings
        string (str): the string to locate

    Returns:
        list
    '''
    assert isinstance(container, (list, tuple))
    ind = []
    for i, s in enumerate(container):
        if string == s:
            ind.append(i)
    return ind


def get_str_indices_by_iden(container, iden):
    '''Return the indices of identified strings in a list or tuple``container``.

    The strings are identified by ``iden``, either a str, int, or a Iterable of these types.
    If ``iden`` is int or corresponding Iterable, the value greater or equal to the
    length of ``container`` will be ignored.

    Args:
        container (list or tuple): container of strings
        iden (int, str, Iterable): the identifier for string to locate

    Returns:
        list, unique indices of identified strings
    '''
    ret = []
    l = len(container)
    if isinstance(iden, int):
        if iden < l:
            ret.append(iden)
    elif isinstance(iden, str):
        ret.extend(get_str_indices(container, iden))
    elif isinstance(iden, Iterable):
        for i in iden:
            if isinstance(i, int):
                if i < l:
                    ret.append(i)
            elif isinstance(i, str):
                ret.extend(get_str_indices(container, i))
    if ret != []:
        return list(OrderedDict.fromkeys(ret).keys())
    return ret


# ====================== PERFORM CALCULATION ======================
# def common_run_calc_cmd(calc_cmd, fout=None, ferr=None):
#     '''
#     Run the calculation command by threading a subprocess calling calc_cmd
#     '''

#     if fout is None:
#         ofile = sp.PIPE
#     else:
#         ofile = open(fout, 'w')

#     if ferr is None:
#         efile = sp.PIPE
#     else:
#         efile = open(ferr, 'w')

#     p = sp.Popen(calc_cmd, stdout=ofile, stderr=efile, shell=True)
#     p.wait()

#     if not fout is None:
#         ofile.close()
#     if not ferr is None:
#         efile.close()


def conv_estimate_number(s):
    '''Convert a string representing a number with error to a float number.

    Literally, string like '3.87(6)' will be converted to 3.876.
    For now, estimate error in the parenthese is reserved.

    Args:
        s (str): number string

    Retuns:
        float
    '''
    return float(re.sub(r"[\(\)]", '', s))


def get_latt_vecs_from_latt_consts(a, b, c, alpha=90, beta=90, gamma=90):
    '''Convert lattice constants to lattice vectors in right-hand system

    Currently support orthormrhobic lattice only!!!

    Args:
        a, b, c (float): length of lattice vectors
        alpha, beta, gamma (float): angles between lattice vectors in degree.
            90 used as default.
    '''
    a = abs(a)
    b = abs(b)
    c = abs(c)
    return [[a, 0, 0], [0, b, 0], [0, 0, c]]


class Cif:
    '''Class to read CIF files and initialize atomic data by PyCIFRW

    Args:
        pathCif (str): the path to cif file
    '''

    def __init__(self, pathCif):
        if not os.path.isfile(pathCif):
            raise FileNotFoundError(pathCif)
        # data block
        self.__blk = CifFile.ReadCif(pathCif, scantype='flex').first_block()
        self.__init_inequiv()
        self.__init_symmetry_operations()
        self._latt = None
        self._atoms = None
        self._pos = None
        self._ref = None

    def __init_inequiv(self):
        '''initialize the positions and symbols of all inequivalent atoms'''
        posInequiv = []
        atomsInequiv = []
        natomsPerInequiv = []
        for l in self.__blk.GetLoop('_atom_site_fract_x'):
            posOne = []
            for a in ["_atom_site_fract_x",
                      "_atom_site_fract_y",
                      "_atom_site_fract_z"]:
                p = conv_estimate_number(l.__getattribute__(a))
                posOne.append(p)
            posInequiv.append(posOne)
            natomsPerInequiv.append(int(l._atom_site_symmetry_multiplicity))
            # remove chemical valence
            atomsInequiv.append(
                re.sub(r"[\d]+[+-]?", '', l._atom_site_type_symbol))
        self._posInequiv = posInequiv
        self._atomsInequiv = atomsInequiv
        self._natomsPerInequiv = natomsPerInequiv

    def __init_symmetry_operations(self):
        '''get all symmetry operations
        '''
        self.operations = {}
        rots = []
        trans = []
        for l in self.__blk.GetLoop('_symmetry_equiv_pos_site_id'):
            r, t = Cif.decode_equiv_pos_string(l._symmetry_equiv_pos_as_xyz)
            rots.append(r)
            trans.append(t)
        self.operations["rotations"] = tuple(rots)
        self.operations["translations"] = tuple(trans)

    def get_lattice_vectors(self):
        '''initialize the lattice vectors from cif

        Returns
            list, shape (3,3)
        '''
        if self._latt is None:
            latta, lattb, lattc = \
                tuple(map(lambda x: conv_estimate_number(self.__blk.GetItemValue(x)),
                        ["_cell_length_a", "_cell_length_b", "_cell_length_c"]))
            angles = []
            for a in ["_cell_angle_alpha",
                    "_cell_angle_beta",
                    "_cell_angle_gamma"]:
                angles.append(conv_estimate_number(self.__blk.GetItemValue(a)))
            self._latt = get_latt_vecs_from_latt_consts(latta, lattb, lattc, *angles)
        return self._latt

    def get_all_atoms(self):
        '''return the symbols and positions of all atoms
        by performing symmetry operations on all inequivalent atoms

        Returns:
            two list, symbols and positions of all atoms, 
            shape (n,) and (n,3) with n the total number of atoms
        '''
        if self._atoms is None or self._pos is None:
            pos = []
            atoms = []
            for r, t in zip(self.operations["rotations"], self.operations["translations"]):
                for i, p in enumerate(self._posInequiv):
                    a = np.add(np.dot(r, p), t)
                    # move to the lattice at origin
                    a = np.subtract(a, np.floor(a))
                    try:
                        for pPrev in pos:
                            if np.allclose(pPrev, a):
                                raise ValueError
                    except ValueError:
                        continue
                    else:
                        atoms.append(self._atomsInequiv[i])
                        pos.append(a)
            # consistency check
            if sum(self._natomsPerInequiv) != len(atoms):
                raise IOError(
                    "inconsistent number of atoms and entries after symmetry operations")
            self._atoms = atoms
            self._pos = pos
        return self._atoms, self._pos
    
    def get_reference_str(self):
        '''Get the reference string

        Returns:
            str
        '''
        if self._ref is None:
            refTitle = ''.join(self.__blk.GetItemValue("_publ_section_title").split('\n'))
            self._ref = refTitle
        return self._ref

    @staticmethod
    def decode_equiv_pos_string(s):
        '''Convert a string representing symmetry operation in CIF file
        to a rotation matrix R and a translation vector t

        The relation between original and transformed fractional coordinate, x and x',
        is

        x' = Rx + t

        Obviously, x, x' and t are treated as a column vector

        Args:
            s (str): a symmetry operation string found in 
                _symmetry_equiv_pos_as_xyz item of a CIF file.

        Returns:
            two lists, shape (3,3) and (3,)
        '''
        trans = [0, 0, 0]
        rot = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        items = [x.strip() for x in s.split(',')]
        if len(items) != 3:
            raise ValueError("s does not seem to be a symmetry operation string")
        for i in items:
            if len(i) == 0:
                raise ValueError(
                    "s does not seem to be a symmetry operation string")

        d = {'x': 0, 'y': 1, 'z': 2}
        for i in range(3):
            stList = items[i].split('+')
            for st in stList:
                # loop in case that string like '-x-y' appears
                while True:
                    sign = 1
                    try:
                        if st.startswith('-'):
                            sign = -1
                            st = st[1:]
                        if st[0] in d:
                            rot[i][d[st[0]]] = sign
                            st = st[1:]
                        else:
                            # confront number
                            break
                    except IndexError:
                        # end of line
                        break
                if len(st) == 0:
                    continue
                else:
                    # deal with fractional number x/y
                    trans[i] = float(st[0]) / float(st[-1])
        return rot, trans
