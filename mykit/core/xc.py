# coding = utf8
'''
'''
import os

from mykit.core._control import (build_tag_map_obj, extract_from_tagdict,
                                 parse_to_tagdict, prog_mapper, tags_mapping)
from mykit.core.log import Verbose


class XCError(Exception):
    pass


class xc_control(Verbose, prog_mapper):
    '''base class that controls exchange-correlation setup
    '''

    # Read tag mapping object from metadata
    _meta = os.path.join(os.path.dirname(__file__), 'metadata', 'xcmap.json')
    _tagMaps = build_tag_map_obj(_meta, "mykit", "json")
    _xcTagMaps = _tagMaps
    _xcValMaps = {}

    def __init__(self, progName, **xctags):
        self._xcTags = {}
        self._parse_xctags(progName, **xctags)

    def parse_tags(self, progName, **xctags):
        '''parse mykit and program-specific tags related to XC

        Note:
            if a program-specific tag and its mykit equivalent
        exists, the program-specific tag value is preferred.

        Args:
            progName : the name of program which the tags should belong to.
        '''
        self._parse_xctags(progName, **xctags)

    def _parse_xctags(self, progName, **xctags):
        if len(xctags) == 0:
            return
        # self.print_log(" In _parse_xctags. Parsing: ", xctags, depth=1, level=3)
        parse_to_tagdict(self._xcTags, self._xcTagMaps, progName, **xctags)
        # self.print_log("End _parse_xctags, now xcTags: ", self._xcTags, depth=1, level=3)
        
    def delete_tags(self, progName, *tags):
        self._pop_xctags(progName, *tags)

    def pop_tags(self, progName, *tags):
        return self._pop_xctags(progName, *tags)

    def _pop_xctags(self, progName, *tags):
        _vals = self._xctag_vals(progName, *tags, delete=True)
        return _vals

    def _get_one_mykit_tag(self, pwTagName):
        '''Get the value of one mykit tag (xc). Conform prog_mapper ABC
        '''
        return self._get_one_xctag(pwTagName)

    def _get_one_xctag(self, xcTagName):
        return self._xcTags.get(xcTagName, None)

    def tag_vals(self, progName, *tags):
        '''find values of xc tags from program-specific tags of progName
        
        The tags name depends on the program, i.e. ``progName``.

        Note:
            Tags that belong to xcTagMaps will also get their value instead of None,
            even when progName is not assigned, i.e. use "mykit"

        Args:
            tags (str): names of tags to request values
            progName (str): the program of which the tags belong to.

        Returns:
            list containing all tag values.
        '''
        return self._xctag_vals(progName, *tags)

    def _xctag_vals(self, progName, *tags, delete=False):
        if len(tags) == 0:
            return []
        # self.print_log("In _xctag_vals, search {} tags of {}: ".format(len(tags), progName), tags, level=3, depth=1)
        vals = extract_from_tagdict(xc_control, self._xcTags, progName, *tags, delete=delete)
        # self.print_log("Found values:", vals, level=3, depth=3)
        return vals

    @classmethod
    def map_tags(cls, *tags, progFrom="mykit", progTo="mykit", getAll=False):
        _pF = progFrom.lower()
        _pT = progTo.lower()
        return tags_mapping(cls._xcTagMaps, _pF, _pT, *tags, getAll=getAll)

    @property
    def xcTags(self):
        return self._xcTags
