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

    # def __initialize_pwtags(self, **kwargs):
    #     '''initialize the fundamental tags of pw class

    #     This will purge out all non-available tags
    #     '''
    #     for _pwt in kwargs:
    #         self.__parse_pwtags(_pwt, kwargs[_pwt])

    def parse_tags(self, **keyval):
        self.__parse_pwtags(**keyval)

    def __parse_pwtags(self, **keyval):
        if len(keyval) == 0:
            return
        for _k, _v in keyval.items():
            if _k in self.__pwTagMaps:
                self.__pwTags.update({_k:_v})
        # if pwTagName in self.__pwTagMaps:
            # self.__pwTags.update({pwTagName: value})

    def __get_one_pwtag(self, pwTagName):
        return self.__pwTags.get(pwTagName, None)

    def tag_vals(self, *tags):
        return self.__tag_vals(*tags, progName="n a")

    def __tag_vals(self, *tags, progName="n a"):
        '''find values of plane_wave tags from program-specific tags of progName
        
        The tags name depends on the program, i.e. ``progName``.

        Returns:
            dict, with key-value pair as : ``program-specific tag: tag value``
        '''
        if len(tags) == 0:
            return None
        _vals = []
        # _pwtags = plane_wave.map2pwtags(*tags, progFrom=progName)
        _pwtags = plane_wave.map_tags(*tags, progFrom=progName, progTo="n a")
        self.print_log("Extracting plane_wave tags: ", _pwtags, level=3, depth=1)
        if len(tags) == 1:
            return self.__get_one_pwtag(_pwtags)
        return list(map(self.__get_one_pwtag, _pwtags))

    @property
    def pwTagvals(self):
        return self.__pwTags

    @classmethod
    def map_tags(cls, *tags, progFrom="n a", progTo="n a"):
        '''classmethod to map tags of progFrom to corresponding tags of progTo

        Args:
            tags (str) : the tag name of progFrom to map from
            progFrom (str)
            progTo (str)
        '''
        # ! Problem in get value of single tag 
        _pF = progFrom.lower()
        _pT = progTo.lower()
        _d = {}
        for _pwt, _map in cls.__pwTagMaps.items():
            # _d.update({cls.__pwTagMaps.get(_pwt, {}).get(_pF, None): \
            #            cls.__pwTagMaps.get(_pwt, {}).get(_pT, None)})
            _d.update({_map.get(_pF, None): _map.get(_pT, None)})
        # ensure tags not implemented will be mapped to None
        _d.update({None:None})
        if len(tags) == 0:
            return tuple(_d.values())
        if len(tags) == 1:
            return _d.get(tags[0], None)
        return tuple(_d.get(_t, None) for _t in tags)

    @classmethod
    def map_from_pwtags(cls, *pwTagNames, progTo="n a"):
        
        '''Class method to map ``plane_wave`` tags to program-specific tags 

        If the pwTagName is not found, None will be return

        Args:
            pwTagNames (str) : the plane_wave tag name to map from
            progName (str) : the name of the programe to which the plane_wave tags will be mapped

        Return:
            a dict value, or a list of dict value corresponding to tagNames
            If no tagNames is specified, all available program-specifi tags will be returned.
        '''
        return cls.map_tags(*pwTagNames, progTo=progTo)

    @classmethod
    def map2pwtags(cls, *tagNames, progFrom="n a"):
        '''Class method to map program-specific tags to ``plane_wave`` tags

        If the tag name is not found, None will be return

        Args:
            tagNames (str) : the name of program-specific tag
            progName (str) : the name of the programe

        Return:
            a dict value, or a list of dict value corresponding to tagNames
            If no tagNames is specified, all available plane_wave tags will be returned.
        '''
        return cls.map_tags(*tagNames, progFrom=progFrom)
