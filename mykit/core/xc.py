# coding = utf8
'''
'''
from mykit.core.log import verbose
from mykit.core._control import control_map

class xcError(Exception):
    pass


class xc_control(verbose, control_map):
    '''base class that controls exchange-correlation setup
    '''

    __xcTagMaps = {
                    "gga" : {"n a":"gga", "vasp":"GGA"},
                    "metagga" : {"n a":"metagga", "vasp":"METAGGA"}
                  }
    __xcValMaps = {}
    __xcTags = {}

    def __init__(self, progName, **xctags):
        self.__parse_xctags(progName, **xctags)

    def parse_tags(self, progName, **xctags):
        self.__parse_xctags(progName, **xctags)

    def __parse_xctags(self, progName, **xctags):
        if len(xctags) == 0:
            return
        self.print_log(" In __parse_xctags. Parsing: ", xctags, depth=1, level=3)
        for _origTag, _v in xctags.items():
            if _origTag == None:
                continue
            elif _origTag in self.__xcTagMaps:
                self.__xcTags.update({_origTag:_v})
            else:
                for _xct, _xcmap in self.__xcTagMaps.items():
                    if _origTag == _xcmap.get(progName, None):
                        self.__xcTags.update({_xct: _v})
                        break
        self.print_log("End __parse_xctags, now xcTags: ", self.__xcTags, depth=1, level=3)
        
    def delete_tags(self, progName, *tags):
        self.__pop_xctags(progName, *tags)

    def pop_tags(self, progName, *tags):
        return self.__pop_xctags(progName, *tags)

    def __pop_xctags(self, progName, *tags):
        _vals = self.__xctag_vals(progName, *tags, delete=True)
        return _vals

    def __get_one_xctag(self, xcTagName):
        return self.__xcTags.get(xcTagName, None)

    def tag_vals(self, progName, *tags):
        '''find values of xc tags from program-specific tags of progName
        
        The tags name depends on the program, i.e. ``progName``.

        Note:
            Tag in ``tags`` that belong to xcTags will get its value instead of None,
            even when progName is not assigned ("n a")

        Args:
            tags (str): names of tags to request values
            progName (str): the program of which the tags belong to.

        Returns:
            list containing all tag values.
        '''
        return self.__xctag_vals(progName, *tags)

    def __xctag_vals(self, progName, *tags, delete=False):
        if len(tags) == 0:
            return []
        self.print_log("In __xctag_vals, search {} tags of {}: ".format(len(tags), progName), tags, level=3, depth=1)
        _xctags = xc_control.map2xctags(*tags, progFrom=progName)
        # self.print_log("_xctags:", _xctags, level=3, depth=2)
        _vals = list(map(self.__get_one_xctag, _xctags))
        # self.print_log("_values:", _vals, level=3, depth=2)
        if progName != "n a":
            for _i, _v in enumerate(_vals):
                if _v == None:
                    if tags[_i] in self.__xcTags:
                        _vals[_i] = self.__xcTags[tags[_i]]
        if delete:
            for _i, _v in enumerate(_xctags):
                if _v in self.__xcTags:
                    del self.__xcTags[_v]
                if tags[_i] in self.__xcTags:
                    del self.__xcTags[tags[_i]]
        self.print_log("Found values:", _vals, level=3, depth=3)
        return _vals

    @classmethod
    def map_tags_in_xc(cls, *tags, progFrom="n a", progTo="n a", getAll=False):
        # if len(tags) == 0:
        #     return tuple()
        _pF = progFrom.lower()
        _pT = progTo.lower()
        return control_map._tags_mapping(cls.__xcTagMaps, _pF, _pT, *tags, getAll=getAll)

    
    @classmethod
    def map2xctags(cls, *tags, progFrom="n a", getAll=False):
        return cls.map_tags_in_xc(*tags, progFrom=progFrom, progTo="n a", getAll=getAll)

    @classmethod
    def map_from_xctags(cls, *tags, progTo="n a", getAll=False):
        return cls.map_tags_in_xc(*tags, progTo=progTo, progFrom="n a", getAll=getAll)

    @property
    def xcTags(self):
        return self.__xcTags