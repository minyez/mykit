# -*- coding: utf-8 -*-
'''define classes and functions related to plane wave basis setup
'''
import os

from mykit.core._control import (build_tag_map_obj, extract_from_tagdict,
                                 parse_to_tagdict, prog_mapper, tags_mapping)
from mykit.core.log import verbose


class PlanewaveError(Exception):
    pass


class planewave_control(verbose, prog_mapper):
    '''the base class that manage parameters of plane-wave basis.
    '''

    # ! Plans of tag mapping of planewave go below. May move to document someday
    # _tagMaps = {
    #             "encutPw":{"mykit":"encutPw", "vasp": "ENCUT"},
    #             # "encutPwGw": {"mykit":"encutPwGw", "vasp": "ENCUTGW"},
    #             # "restartWave":{"mykit": "restartWave", "vasp": "ISTART"},
    #             # "restartCharg":{"mykit": "restartCharg", "vasp": "ICHARG"},
    #             # "nscf": {"mykit": "nscf", "vasp": "NELM"},
    #             # "ediff": {"mykit": "ediff", "vasp": "EDIFF"},
    #             # "scfAlgo": {"mykit": "scfAlgo", "vasp": "ALGO"},
    #             # "globalPrec": {"mykit": "globalPrec", "vasp": "PREC"},
    #             # "ifWriteWave": {"mykit": "ifWriteWave", "vasp": "LWAVE"},
    #             # "ifWriteCharg": {"mykit": "ifWriteCharg", "vasp": "LCHARG"},
    #            }
    # Read tag mapping object from metadata
    _meta = os.path.join(os.path.dirname(__file__), 'metadata', 'planewavemap.json')
    _tagMaps = build_tag_map_obj(_meta, "mykit", "json")
    _pwTagMaps = _tagMaps
    _pwValMaps = {}

    def __init__(self, progName, **pwargs):
        self._pwTags = {}
        # if len(pwargs) > 0:
        self._parse_pwtags(progName, **pwargs)
        # else:
        #     pass

    def parse_tags(self, progName, **pwtags):
        '''parse mykit and program-specific tags relatex to plane waves.

        Note:
            if a program-specific tag and its mykit equivalent
        exists, the program-specific tag value is preferred.

        Args:
            progName : the name of program which the tags should belong to.
        '''
        self._parse_pwtags(progName, **pwtags)

    def _parse_pwtags(self, progName, **pwtags):
        if len(pwtags) == 0:
            return
        # self.print_log(" In _parse_pwtags. Parsing {} Tags ".format(progName), pwtags, depth=1, level=3)
        # self.print_log("                    pwTags before:", self._pwTags, depth=1, level=3)
        parse_to_tagdict(self._pwTags, self._pwTagMaps, progName, **pwtags)
        # self.print_log("End _parse_pwtags. pwTags after:", self._pwTags, depth=1, level=3)

    def delete_tags(self, progName, *tags):
        self._pop_pwtags(progName, *tags)

    def pop_tags(self, progName, *tags):
        return self._pop_pwtags(progName, *tags)

    def _pop_pwtags(self, progName, *tags):
        _vals = self._pwtag_vals(progName, *tags, delete=True)
        return _vals

    def _get_one_mykit_tag(self, pwTagName):
        return self._get_one_pwtag(pwTagName)

    def _get_one_pwtag(self, pwTagName):
        return self._pwTags.get(pwTagName, None)

    def tag_vals(self, progName, *tags):
        '''find values of plane_wave tags from program-specific tags of progName
        
        The tags name depends on the program, i.e. ``progName``.

        Note:
            Tag in ``tags`` that belong to pwTags will get its value instead of None,
            even when progName is not assigned ("mykit")

        Args:
            tags (str): names of tags to request values
            progName (str): the program of which the tags belong to.

        Returns:
            list containing all tag values.
        '''
        return self._pwtag_vals(progName, *tags)

    def _pwtag_vals(self, progName, *tags, delete=False):
        if len(tags) == 0:
            return []
        # self.print_log("In _pwtag_vals, search {} tags of {}: ".format(len(tags), progName), tags, level=3, depth=1)
        vals = extract_from_tagdict(planewave_control, self._pwTags, progName, *tags, delete=delete)
        # self.print_log("Found values:", vals, level=3, depth=3)
        return vals

    @property
    def pwTags(self):
        return self._pwTags

    @classmethod
    def map_tags(cls, *tags, progFrom="mykit", progTo="mykit", getAll=False):
        '''classmethod to map tags of progFrom to corresponding tags of progTo

        Args:
            tags (str) : the tag name of progFrom to map from
            progFrom (str)
            progTo (str)
            getAll (bool)
        '''
        _pF = progFrom.lower()
        _pT = progTo.lower()
        return tags_mapping(cls._pwTagMaps, _pF, _pT, *tags, getAll=getAll)
