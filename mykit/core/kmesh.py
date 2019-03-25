# -*- coding: utf-8 -*-
'''Module that defines classes and functions for Brillouin zone sampling
'''
import os
import re
from copy import deepcopy

import numpy as np

from mykit.core._control import (build_tag_map_obj, extract_from_tagdict,
                                 parse_to_tagdict, prog_mapper, tags_mapping)
from mykit.core.log import Verbose
from mykit.core.numeric import Prec

# from mykit.core.utils import if_vec_same_direction

# Allowed pattern for kpoint symbols and kpath string
KSYM_PATTERN = r'[A-Z]{1,2}'
KPATH_PATTERN = r'^('+ KSYM_PATTERN + r'-)+' + KSYM_PATTERN + r'$'

class KmeshError(Exception):
    pass


class kmesh_control(Verbose, prog_mapper):

    _meta = os.path.join(os.path.dirname(__file__), 'metadata', 'kmeshmap.json')
    _tagMaps = build_tag_map_obj(_meta, "mykit", "json")
    _kmeshTagMaps = _tagMaps
    _kmeshValMaps = {}

    def __init__(self, progName, **kmargs):
        self._kmeshTags = {}
        self._parse_kmeshtags(progName, **kmargs)
    
    def parse_tags(self, progName, **kmtags):
        '''
        '''
        self._parse_kmeshtags(progName, **kmtags)
    
    def _parse_kmeshtags(self, progName, **kmtags):
        if len(kmtags) == 0:
            return
        parse_to_tagdict(self._kmeshTags, self._kmeshTagMaps, progName, **kmtags)
    
    def delete_tags(self, progName, *tags):
        self._pop_kmeshtags(progName, *tags)
    
    def pop_tags(self, progName, *tags):
        return self._pop_kmeshtags(progName, *tags)

    def _pop_kmeshtags(self, progName, *tags):
        vals = self._kmeshtag_vals(progName, *tags, delete=True)
        return vals

    def _get_one_mykit_tag(self, kmTagName):
        return self._get_one_kmeshtag(kmTagName)
    
    def _get_one_kmeshtag(self, kmTagName):
        return self._kmeshTags.get(kmTagName, None)
    
    def tag_vals(self, progName, *tags):
        return self._kmeshtag_vals(progName, *tags)
    
    def _kmeshtag_vals(self, progName, *tags, delete=False):
        if len(tags) == 0:
            return []
        vals = extract_from_tagdict(kmesh_control, self._kmeshTags, progName, *tags, delete=delete)
        return vals
    
    @property
    def kmeshTags(self):
        return self._kmeshTags

    @property
    def kmode(self):
        return self._kmeshTags.get("kmode")

    @property
    def kdiv(self):
        return self._kmeshTags.get("div")

    @classmethod
    def map_tags(cls, *tags, progFrom="mykit", progTo="mykit", getAll=False):
        '''
        '''
        _pF = progFrom.lower()
        _pT = progTo.lower()
        return tags_mapping(cls._kmeshTagMaps, _pF, _pF, *tags, getAll=getAll)


def kpath_decoder(kpath):
    '''Decode string consisting kpoints symbol to a list of strings.

    Those with even and odd indices (from 0) are the starting and 
    ending point in reciprocal space, repectively

    The path can have more than one continuous path in reciprocal space,
    separated by space.
    However, each continuous path should match the ``KPATH_PATTERN``, 
    otherwise `KmeshError` will be raised.

    Args:
        kpath (str): the string containing kpoints symbol and 
        representing a trajectory in reciprocal space

    Examples:
    >>> kpath_decoder("A-B-C D-E")
    ["A", "B", "B", "C", "D", "E"]
    >>> kpath_decoder("GM-W-X-L-GM-X")
    ["GM", "W", "W", "X", "X", "L", "L", "GM", "GM", "X"]
    '''
    try:
        _klines = kpath.split()
    except (AttributeError, SyntaxError):
        raise KmeshError("Input kpath should be string: {}".format(kpath))

    # the pattern of each path segment
    linePat = re.compile(KPATH_PATTERN)

    ksegs = []
    for kline in _klines:
        if not re.match(linePat, kline):
            raise KmeshError("Invalid kpath line string: {}".format(kline))
        symbols = kline.split('-')
        nSyms = len(symbols)
        for i in range(nSyms-1):
            # ksegs.append('{}-{}'.format(symbols[i], symbols[i+1]))
            if symbols[i] == symbols[i+1]:
                raise KmeshError("kpath with zero length: {}-{}".format(symbols[i], symbols[i+1]))
            ksegs.extend([symbols[i], symbols[i+1]])
    return ksegs


def kpath_encoder(ksyms):
    '''Encode a list/tuple of strings to a complete kpath string.

    Args:
        ksyms (list or tuple): container of kpath symbols, must have an even length
    '''
    try:
        assert isinstance(ksyms, (list, tuple))
    except AssertionError:
        raise KmeshError("require list or tuple, received {}".format(type(ksyms)))
    try:
        assert len(ksyms)%2 == 0
    except AssertionError:
        raise KmeshError("require even length, received {}".format(len(ksyms)))

    nLineSeg = int(len(ksyms)/2)
    kpath = ''
    symPat = re.compile(r'^' + KSYM_PATTERN + r'$')
    lastSym = ''

    for _i in range(nLineSeg):
        st = ksyms[2*_i]
        ed = ksyms[2*_i+1]
        if st == ed:
            raise KmeshError("kpath with zero length: {}-{}".format(st, ed))
        seg = (st, ed)
        for ksym in seg:
            if not re.match(symPat, ksym):
                raise KmeshError("Invalid kpoint symbol: {}".format(ksym))
        if _i == 0:
            kpath += '-'.join(seg)
        else:
            if st == lastSym:
                kpath += '-' + ed
            else:
                kpath += ' ' + '-'.join(seg)
        lastSym = ed
    return kpath


def _check_valid_ksym_coord_pair(ksym, coord):
    if not re.match(r"^" + KSYM_PATTERN + r"$", ksym):
        raise KeyError("Invalid kpoint symbol: {}".format(ksym))
    try:
        shape = np.shape(coord)
    except ValueError:
        raise ValueError("Invalid kpoint coordinate for symbol {}".format(ksym))
    else:
        if shape != (3,):
            raise ValueError("Invalid kpoint coordinate for symbol {}".format(ksym))

def _check_valid_kpath_dict(kpathDict):
    try:
        assert isinstance(kpathDict, dict)
    except AssertionError:
        raise TypeError("kpath must be dictionary.")
    try:
        assert set(["symbols", "coordinates"]) == set(kpathDict.keys())
    except AssertionError:
        raise KeyError("\"symbols\", \"coordinates\" keys not found. Please check")
    for (ksym, coord) in zip(kpathDict["symbols"],kpathDict["coordinates"]):
        try:
            _check_valid_ksym_coord_pair(ksym, coord)
        except (KeyError, ValueError) as _err:
            raise _err

def check_kvecs_form_kpath(kvec):
    '''Check if the kpoint vectors form several line segments in the reciprocal space

    Usually, the number of kpoints on one line segments is no less than 3.

    Args:
        kvec (array-like): the kpoint vectors to analysis, shape, (n,3) 

    Returns:
        list, with tuple as members. Each tuple has 2 int members,
        the indices of kpoint vectors at the beginning and end of 
        a line segment
    '''
    segs = []
    # check the shape of kvec
    try:
        shape = np.shape(kvec)
        assert len(shape) == 2
        assert shape[1] == 3
    except (TypeError, AssertionError):
        return segs
    nkpt = shape[0]
    if nkpt < 3:
        return segs

    _kvec = np.array(kvec, dtype=Prec._dtype)
    dkvec = _kvec[1:, :] - _kvec[:-1, :]
    # normalize the dkvec, 
    n = np.linalg.norm(dkvec, axis=1)
    for i in range(nkpt-1):
        if np.isclose(n[i], 0):
            dkvec[i,:] = 1000.0
        else:
            dkvec[i,:] = dkvec[i,:]/n[i]
    dp = np.sum(dkvec[:-1,:] * dkvec[1:,:], axis=1)
    st = 0
    ed = 2
    while ed < nkpt:
        if not np.isclose(dp[ed-2], 1):
            if ed - st > 2:
                segs.append((st, ed-1))
            st = ed - 1
        ed = ed + 1
    if ed - st > 2:
        segs.append((st, ed-1))
    return segs

# the mapping from kpoint symbol to LaTeX commands
#pylint: disable=anomalous-backslash-in-string
KSYMBOL_LATEX = {
    "GM": "$\Gamma$",
    "LM": "$\lambda$",
}
