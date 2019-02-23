# -*- coding: utf-8 -*-
'''define classes and functions related to plane wave basis setup
'''
from mykit.core.log import verbose

class planewaveError(Exception):
    pass


class plane_wave(verbose):
    '''the base class that manage parameters of plane-wave basis.
    '''

    __pwTagMaps = {
                    "encutPw":{"n a":"encutPw", "vasp": "ENCUT"},
                    "encutPwGw": {"n a":"encutPwGw", "vasp": "ENCUTGW"},
                    "restartWave":{"n a": "restartWave", "vasp": "ISTART"},
                    "restartCharg":{"n a": "restartCharg", "vasp": "ICHARG"},
                    "nscf": {"n a": "nscf", "vasp": "NELM"},
                    "ediff": {"n a": "ediff", "vasp": "EDIFF"},
                    "scfAlgo": {"n a": "scfAlgo", "vasp": "ALGO"},
                    "globalPrec": {"n a": "globalPrec", "vasp": "PREC"},
                    "ifWriteWave": {"n a": "ifWriteWave", "vasp": "LWAVE"},
                    "ifWriteCharg": {"n a": "ifWriteCharg", "vasp": "LCHARG"},
                   }
    __pwTags = {}

    def __init__(self, **pwargs):
        if len(pwargs) > 0:
            self.__parse_pwtags(**pwargs)
        else:
            pass

    def parse_tags(self, **keyval):
        self.__parse_pwtags(**keyval)

    def __parse_pwtags(self, **keyval):
        if len(keyval) == 0:
            return
        for _k, _v in keyval.items():
            if _k in self.__pwTagMaps:
                self.__pwTags.update({_k:_v})

    def __get_one_pwtag(self, pwTagName):
        return self.__pwTags.get(pwTagName, None)

    def tag_vals(self, *tags, progName="n a"):
        return self.__pwtag_vals(*tags, progName=progName)

    def __pwtag_vals(self, *tags, progName="n a"):
        '''find values of plane_wave tags from program-specific tags of progName
        
        The tags name depends on the program, i.e. ``progName``.

        Returns:
            list, if tags is specified, otherwise None
        '''
        if len(tags) == 0:
            return None
        # _vals = []
        _pwtags = plane_wave.map2pwtags(*tags, progFrom=progName)
        # self.print_log("Extracting plane_wave tags: ", _pwtags, level=3, depth=1)
        # ? get value from plane_wave tag, even progName is not "n a"
        _vals = list(map(self.__get_one_pwtag, _pwtags))
        if progName != "n a":
            for _i, _v in enumerate(_vals):
                if _v == None:
                    if tags[_i] in self.__pwTags:
                        _vals[_i] = self.__pwTags[tags[_i]]
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
        _d = {}
        for _pwt, _map in cls.__pwTagMaps.items():
            _d.update({_map.get(_pF, None): _map.get(_pT, None)})
        # ensure tags not implemented will be mapped to None
        _d.update({None:None})
        if getAll:
            _d.pop(None)
            return tuple(_d.values())
        if len(tags) == 0:
            return tuple()
        if len(tags) == 1:
            return (_d.get(tags[0], None),)
        return tuple(_d.get(_t, None) for _t in tags)

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
