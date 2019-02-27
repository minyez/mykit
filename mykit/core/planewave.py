# -*- coding: utf-8 -*-
'''define classes and functions related to plane wave basis setup
'''
from mykit.core.log import verbose
from mykit.core._control import tags_mapping, prog_mapper

class PlanewaveError(Exception):
    pass


class planewave_control(verbose, prog_mapper):
    '''the base class that manage parameters of plane-wave basis.

    parse_tags method for all base class need to specify a progName argument
    '''

    _pwTagMaps = {
                    "encutPw":{"mykit":"encutPw", "vasp": "ENCUT"},
                    # "encutPwGw": {"mykit":"encutPwGw", "vasp": "ENCUTGW"},
                    # "restartWave":{"mykit": "restartWave", "vasp": "ISTART"},
                    # "restartCharg":{"mykit": "restartCharg", "vasp": "ICHARG"},
                    # "nscf": {"mykit": "nscf", "vasp": "NELM"},
                    # "ediff": {"mykit": "ediff", "vasp": "EDIFF"},
                    # "scfAlgo": {"mykit": "scfAlgo", "vasp": "ALGO"},
                    # "globalPrec": {"mykit": "globalPrec", "vasp": "PREC"},
                    # "ifWriteWave": {"mykit": "ifWriteWave", "vasp": "LWAVE"},
                    # "ifWriteCharg": {"mykit": "ifWriteCharg", "vasp": "LCHARG"},
                   }
    _pwValMaps = {}
    _pwTags = {}

    def __init__(self, progName, **pwargs):
        # if len(pwargs) > 0:
        self._parse_pwtags(progName, **pwargs)
        # else:
        #     pass

    def parse_tags(self, progName, **pwtags):
        '''parse plane_wave and program-specific tags to pwTags.

        Note:
            if a program-specific tag and its plane_wave correspondent
        exists, the program-specific tag value is preferred.

        Args:
            progName : the name of program that is required.
        '''
        self._parse_pwtags(progName, **pwtags)

    def _parse_pwtags(self, progName, **pwtags):
        if len(pwtags) == 0:
            return
        self.print_log(" In _parse_pwtags. Parsing {} Tags ".format(progName), pwtags, depth=1, level=3)
        self.print_log("                    pwTags before:", self._pwTags, depth=1, level=3)
        for _origTag, _v in pwtags.items():
            if _origTag == None:
                continue
            # check plane_wave tags
            elif _origTag in self._pwTagMaps.keys():
                self._pwTags.update({_origTag:_v})
            else:
            # check program-specific tags
            # ? Can be optimized
                for _pwt, _pwmap in self._pwTagMaps.items():
                    if _origTag == _pwmap.get(progName, None):
                        self._pwTags.update({_pwt: _v})
                        break
        self.print_log("End _parse_pwtags. pwTags after:", self._pwTags, depth=1, level=3)

    def delete_tags(self, progName, *tags):
        self._pop_pwtags(progName, *tags)

    def pop_tags(self, progName, *tags):
        return self._pop_pwtags(progName, *tags)

    def _pop_pwtags(self, progName, *tags):
        _vals = self._pwtag_vals(progName, *tags, delete=True)
        return _vals

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
        self.print_log("In _pwtag_vals, search {} tags of {}: ".format(len(tags), progName), tags, level=3, depth=1)
        _pwtags = planewave_control.map_to_mykit_tags(*tags, progFrom=progName)
        # self.print_log("Extracting plane_wave tags: ", _pwtags, level=3, depth=1)
        # ? get value from plane_wave tag, even progName is not "mykit"
        _vals = list(map(self._get_one_pwtag, _pwtags))
        # No need to get pwtags again, when program is not assigned
        if progName != "mykit":
            for _i, _v in enumerate(_vals):
                if _v == None:
                    if tags[_i] in self._pwTags:
                        _vals[_i] = self._pwTags[tags[_i]]
        if delete:
            for _i, _v in enumerate(_pwtags):
                if _v in self._pwTags:
                    del self._pwTags[_v]
                if tags[_i] in self._pwTags:
                    del self._pwTags[tags[_i]]
        self.print_log("Found values:", _vals, level=3, depth=3)
        return _vals

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
