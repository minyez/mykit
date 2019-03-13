# -*- coding: utf-8 -*-
'''Module that defines classes and functions for Brillouin zone sampling
'''
import os
import re

from mykit.core._control import (build_tag_map_obj, extract_from_tagdict,
                                 parse_to_tagdict, prog_mapper, tags_mapping)
from mykit.core.log import verbose

# from collections.abc import Iterable



class KmeshError(Exception):
    pass


class kmesh_control(verbose, prog_mapper):

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
    def kgrid(self):
        return self._kmeshTags.get("kgrid")

    @classmethod
    def map_tags(cls, *tags, progFrom="mykit", progTo="mykit", getAll=False):
        '''
        '''
        _pF = progFrom.lower()
        _pT = progTo.lower()
        return tags_mapping(cls._kmeshTagMaps, _pF, _pF, *tags, getAll=getAll)


def kpath_decoder(kpath):
    '''Decode string consisting kpoints symbol to a list of strings,
    each of which represents a straight line segment in reciprocal space

    The path can have more than one continuous path in reciprocal space,
    separated by space.
    However, each continuous path should match the pattern 
    `^([A-Z]{1,2}-)+[A-Z]{1,2}$`, i.e. words with 1 or 2 captical
    characters, joint by one hyphen. `KmeshError` will be raised if the
    pattern is not matched.

    Args:
        kpath (str): the string containing kpoints symbol and 
        representing a trajectory in reciprocal space

    Examples:
    >>> kpath_decoder("A-B-C D-E")
    ["A-B", "B-C", "D-E"]
    >>> kpath_decoder("GM-W-X-L-GM-X")
    ["GM-W", "W-X", "X-L", "L-GM", "GM-X"]
    '''
    try:
        _klines = kpath.split()
    except (AttributeError, SyntaxError):
        raise KmeshError("Input kpath should be string: {}".format(kpath))

    # the pattern of each path segment
    linePat = re.compile(r'^([A-Z]{1,2}-)+[A-Z]{1,2}$')

    ksegs = []
    for kline in _klines:
        if not re.match(linePat, kline):
            raise KmeshError("Invalid kpath line string: {}".format(kline))
        symbols = kline.split('-')
        nSyms = len(symbols)
        for i in range(nSyms-1):
            ksegs.append('{}-{}'.format(symbols[i], symbols[i+1]))
    return ksegs


def kpath_encoder(ksegs):
    '''Encode a list of kpath segment strings to a complete kpath string.

    Args:
        ksegs (list or tuple): container of kpath segment strings
    '''
    try:
        assert isinstance(ksegs, (list, tuple))
    except:
        raise KmeshError("require list or tuple, received {}".format(type(ksegs)))

    kpath = ''
    segPat = re.compile(r'^[A-Z]{1,2}-[A-Z]{1,2}$')
    lastSym = ''

    for i, kseg in enumerate(ksegs):
        if not re.match(segPat, kseg):
            raise KmeshError("Invalid kpath segment string: {}".format(kseg))
        if i == 0:
            kpath += kseg
            lastSym = kseg.split('-')[-1]
            continue
        syms = kseg.split('-')
        if syms[0] == lastSym:
            kpath += '-' + syms[1]
        else:
            kpath += ' ' + kseg
        lastSym = kseg.split('-')[-1]
    return kpath
