# -*- coding: utf-8 -*-
'''Module that defines classes and functions for Brillouin zone sampling
'''
import os

from mykit.core._control import (build_tag_map_obj, extract_from_tagdict,
                                 parse_to_tagdict, prog_mapper, tags_mapping)
from mykit.core.log import verbose


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
