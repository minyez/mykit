# -*- coding: utf-8 -*-
'''define classes and functions related to plane wave basis setup
'''

class plane_wave:
    '''the class that defines variables that control the plane wave basis

    This is the base class for classes that control input parameters for calculations
    in plane-wave software.
    '''

    __progName = "n a" # not assigned
    __encutPw = None
    __restart = 0
    __encutPwGw = None
    __mapProgTag = {
                    "encutPw":{"n a":"encutPw", "vasp": "ENCUT"},
                    "encutPwGw": {"n a":"encutPwGw", "vasp": "ENCUTGW"},
                    "restart":{"n a": "restart", "vasp": "ISTART"},
                   }
    __pwTags = {}

    def __init__(self, program, **kwargs):
        __progName = program.lower()
        if len(kwargs) > 0:
            self.__pwTags = kwargs
        self.__initialize_pwtags()
    
    def __initialize_pwtags(self):
        '''initialize the fundamental tags of pw class

        This will purge out all non-available tags
        '''
        for _pwt, _v in enumerate(self.__pwTags):
            if _pwt in self.__mapProgTag:
                self._parse_one_pwtag(_pwt, _v)
            else:
                self.__pwTags.pop(_pwt)

    def _parse_one_pwtag(self, pwTagName, value):
        if pwTagName == "encutPw":
            self.__encutPw = value
        elif pwTagName == "encutPwGw":
            self.__encutPwGw = value
        elif pwTagName == "restart":
            self.__restart = value

    def __get_pwtag_value(self, pwTagName):
        return self.__pwTags.get(pwTagName, None)

    # def __back_map_pw_tag(self, tag):
    #     '''back map the tag for particular program to the pw tag
    #     '''
    #     pass

    @property
    def progName(self):
        '''Return the name of the program which the plane_wave class control is currently for
        '''
        return self.__progName

    @property
    def pwTagvals(self):
        return self.__pwTags

    def tagvals(self):
        '''All values of plane_wave tags that have been set. 
        
        The tags name depends on the program, i.e. ``progName``.

        Returns:
            dict, with key-value pair as : ``program-specific tag: tag value``
        '''
        _d = {}
        for _pwt in self.__pwTags:
            _tDict = self.__mapProgTag.get(_pwt, None)
            if _tDict is not None:
                _progTag = _tDict.get(self.__progName, None)
                if _progTag is not None:
                    _d.update({_progTag: self.__pwTags[_pwt]})
        return _d

    # @classmethod
    # def __map_one_tag(cls, )

    @classmethod
    def map_tags(cls, *tags, progFrom="n a", progTo="n a"):
        '''classmethod to map tags of progFrom to corresponding tags of progTo

        Args:
            tags (str) : the tag name of progFrom to map from
            progFrom (str)
            progTo (str)
        '''
        _pF = progFrom.lower()
        _pT = progTo.lower()
        _d = {}
        
        for _pwt in cls.__mapProgTag:
            _d.update({cls.__mapProgTag.get(_pwt, {}).get(_pF, None): \
                       cls.__mapProgTag.get(_pwt, {}).get(_pT, None)})
        # ensure tags not implemented will be mapped to None
        _d.update({None:None})
        if len(tags) == 0:
            return tuple(_d.values())
        if len(tags) == 1:
            return _d.get(tags, None)
        return tuple(_d.get(_t, None) for _t in tags)

    @classmethod
    def mapFromPwTags(cls, *pwTagNames, progName="n a"):
        
        '''Class method to map ``plane_wave`` tags to program-specific tags 

        If the pwTagName is not found, None will be return

        Args:
            pwTagNames (str) : the plane_wave tag name to map from
            progName (str) : the name of the programe to which the plane_wave tags will be mapped

        Return:
            a dict value, or a list of dict value corresponding to tagNames
            If no tagNames is specified, all available program-specifi tags will be returned.
        '''
        return cls.map_tags(*pwTagNames, progTo=progName)

    @classmethod
    def map2PwTags(cls, *tagNames, progName="n a"):
        '''Class method to map program-specific tags to ``plane_wave`` tags

        If the tag name is not found, None will be return

        Args:
            tagNames (str) : the name of program-specific tag
            progName (str) : the name of the programe

        Return:
            a dict value, or a list of dict value corresponding to tagNames
            If no tagNames is specified, all available plane_wave tags will be returned.
        '''
        return cls.map_tags(*tagNames, progFrom=progName)
