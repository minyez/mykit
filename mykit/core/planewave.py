# -*- coding: utf-8 -*-
'''define classes and functions related to plane wave basis setup
'''
from mykit.core.log import verbose
from mykit.core._control import control_map

class planewaveError(Exception):
    pass


class plane_wave_control(verbose, control_map):
    '''the base class that manage parameters of plane-wave basis.

    parse_tags method for all base class need to specify a progName argument
    '''

    __pwTagMaps = {
                    "encutPw":{"n a":"encutPw", "vasp": "ENCUT"},
                    "encutPwGw": {"n a":"encutPwGw", "vasp": "ENCUTGW"},
                    "restartWave":{"n a": "restartWave", "vasp": "ISTART"},
                    "restartCharg":{"n a": "restartCharg", "vasp": "ICHARG"},
                    # "nscf": {"n a": "nscf", "vasp": "NELM"},
                    # "ediff": {"n a": "ediff", "vasp": "EDIFF"},
                    # "scfAlgo": {"n a": "scfAlgo", "vasp": "ALGO"},
                    # "globalPrec": {"n a": "globalPrec", "vasp": "PREC"},
                    # "ifWriteWave": {"n a": "ifWriteWave", "vasp": "LWAVE"},
                    # "ifWriteCharg": {"n a": "ifWriteCharg", "vasp": "LCHARG"},
                   }
    __pwValMaps = {}
    __pwTags = {}

    def __init__(self, progName, **pwargs):
        # if len(pwargs) > 0:
        self.__parse_pwtags(progName, **pwargs)
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
        self.__parse_pwtags(progName, **pwtags)

    def delete_tags(self, progName, *tags):
        self.__pop_pwtags(progName, *tags)

    def pop_tags(self, progName, *tags):
        return self.__pop_pwtags(progName, *tags)

    def __pop_pwtags(self, progName, *tags):
        _vals = self.__pwtag_vals(progName, *tags, delete=True)
        return _vals

    def __parse_pwtags(self, progName, **pwtags):
        if len(pwtags) == 0:
            return
        for _origTag, _v in pwtags.items():
            if _origTag == None:
                continue
            # check plane_wave tags
            elif _origTag in self.__pwTagMaps:
                self.__pwTags.update({_origTag:_v})
            else:
            # check program-specific tags
            # ? Can be optimized
                for _pwt, _pwmap in self.__pwTagMaps.items():
                    if _origTag == _pwmap.get(progName, None):
                        self.__pwTags.update({_pwt: _v})
                        break
        # self.print_log("__parse_pwtags: pwTags", self.__pwTags, depth=1, level=3)

    def __get_one_pwtag(self, pwTagName):
        return self.__pwTags.get(pwTagName, None)

    def tag_vals(self, progName, *tags):
        '''find values of plane_wave tags from program-specific tags of progName
        
        The tags name depends on the program, i.e. ``progName``.

        Note:
            Tag in ``tags`` that belong to pwTags will get its value instead of None,
            even when progName is not assigned ("n a")

        Args:
            tags (str): names of tags to request values
            progName (str): the program of which the tags belong to.

        Returns:
            list containing all tag values.
        '''
        return self.__pwtag_vals(progName, *tags)

    def __pwtag_vals(self, progName, *tags, delete=False):
        if len(tags) == 0:
            return []
        # self.print_log("In __pwtag_vals, search {} tags of {}: ".format(len(tags), progName), tags, level=3, depth=1)
        _pwtags = plane_wave_control.map2pwtags(*tags, progFrom=progName)
        # self.print_log("_pwtags:", _pwtags, level=3, depth=2)
        # self.print_log("Extracting plane_wave tags: ", _pwtags, level=3, depth=1)
        # ? get value from plane_wave tag, even progName is not "n a"
        _vals = list(map(self.__get_one_pwtag, _pwtags))
        # No need to get pwtags again, when program is not assigned
        if progName != "n a":
            for _i, _v in enumerate(_vals):
                if _v == None:
                    if tags[_i] in self.__pwTags:
                        _vals[_i] = self.__pwTags[tags[_i]]
        if delete:
            for _i, _v in enumerate(_pwtags):
                if _v in self.__pwTags:
                    del self.__pwTags[_v]
                if tags[_i] in self.__pwTags:
                    del self.__pwTags[tags[_i]]
        return _vals

    @property
    def pwTags(self):
        return self.__pwTags

    @classmethod
    def map_tags_in_pw(cls, *tags, progFrom="n a", progTo="n a", getAll=False):
        '''classmethod to map tags of progFrom to corresponding tags of progTo

        Args:
            tags (str) : the tag name of progFrom to map from
            progFrom (str)
            progTo (str)
        '''
        _pF = progFrom.lower()
        _pT = progTo.lower()
        return control_map._tags_mapping(cls.__pwTagMaps, _pF, _pT, *tags, getAll=getAll)

    @classmethod
    def map_from_pwtags(cls, *pwTagNames, progTo="n a", getAll=False):
        
        '''Class method to map ``plane_wave`` tags to program-specific tags 

        If the pwTagName is not found, None will be return

        Args:
            pwTagNames (str) : the plane_wave tag name to map from
            progName (str) : the name of the programe to which the plane_wave tags will be mapped

        Return:
            a dict value, or a list of dict value corresponding to tagNames
            If no tagNames is specified, all available program-specifi tags will be returned.
        '''
        return cls.map_tags_in_pw(*pwTagNames, progTo=progTo, progFrom="n a", getAll=getAll)

    @classmethod
    def map2pwtags(cls, *tagNames, progFrom="n a", getAll=False):
        '''Class method to map program-specific tags to ``plane_wave`` tags

        If the tag name is not found, None will be return

        Args:
            tagNames (str) : the name of program-specific tag
            progName (str) : the name of the programe

        Return:
            a dict value, or a list of dict value corresponding to tagNames
            If no tagNames is specified, all available plane_wave tags will be returned.
        '''
        return cls.map_tags_in_pw(*tagNames, progFrom=progFrom, progTo="n a", getAll=getAll)
