# coding = utf-8

import os

from mykit.core._control import (build_tag_map_obj, extract_from_tagdict,
                                 parse_to_tagdict, prog_mapper, tags_mapping)
from mykit.core.log import verbose


class IonError(Exception):
    pass


class ion_control(verbose, prog_mapper):

    _meta = os.path.join(os.path.dirname(__file__), 'metadata', 'ionmap.json')
    _tagMaps = build_tag_map_obj(_meta, "mykit", "json")
    _ionTagMaps = _tagMaps
    _ionValMaps = {}

    def __init__(self, progName, **iontags):
        self._ionTags = {}
        self._parse_iontags(progName, **iontags)
    
    def parse_tags(self, progName, **iontags):
        self._parse_iontags(progName, **iontags)

    def pop_tags(self, progName, *tags):
        return self._pop_iontags(progName, *tags)

    def delete_tags(self, progName, *tags):
        self._pop_iontags(progName, *tags)

    def tag_vals(self, progName, *tags):
        return self._iontag_vals(progName, *tags)

    def _parse_iontags(self, progName, **iontags):
        if len(iontags) == 0:
            return
        parse_to_tagdict(self._ionTags, self._ionTagMaps, progName, **iontags)

    def _pop_iontags(self, progName, *tags):
        return self._iontag_vals(progName, *tags, delete=True)

    def _get_one_mykit_tag(self, ionTagName):
        return self._get_one_iontag(ionTagName)

    def _get_one_iontag(self, ionTagName):
        return self._ionTags.get(ionTagName, None)

    def _iontag_vals(self, progName, *tags, delete=False):
        if len(tags) == 0:
            return []
        return extract_from_tagdict(ion_control, self._ionTags, progName, *tags, delete=delete)

    @property
    def ionTags(self):
        return self._ionTags

    @classmethod
    def map_tags(cls, *tags, progFrom="mykit", progTo="mykit", getAll=False):
        return tags_mapping(cls._ionTagMaps, progFrom, progTo, *tags, getAll=getAll)
