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
        self.print_log("__parse_xctags: xcTags", self.__xcTags, depth=1, level=3)
        

    def __get_one_xctag(self, xcTagName):
        return self.__xcTags.get(xcTagName, None)

    def tag_vals(self, *tags, progName="n a"):
        return self.__xctag_vals(*tags, progName=progName)

    def __xctag_vals(self, *tags, progName="n a"):
        if len(tags) == 0:
            return []
        self.print_log("In __xctag_vals, search {} tags of {}: ".format(len(tags), progName), tags, level=3, depth=1)
        _xctags = xc_control.map2xctags(*tags, progFrom=progName)
        self.print_log("_xctags:", _xctags, level=3, depth=2)
        _vals = list(map(self.__get_one_xctag, _xctags))
        self.print_log("_values:", _vals, level=3, depth=2)
        if progName != "n a":
            for _i, _v in enumerate(_vals):
                if _v == None:
                    if tags[_i] in self.__xcTags:
                        _vals[_i] = self.__xcTags[tags[_i]]
        return _vals

    @classmethod
    def map_tags_in_xc(cls, *tags, progFrom="n a", progTo="n a", getAll=False):
        # if len(tags) == 0:
        #     return tuple()
        _pF = progFrom.lower()
        _pT = progTo.lower()
        # _d = {}
        # # cls.print_cm_log("In map_tags_in_xc", level=3, depth=1)
        # for _xct, _map in cls.__xcTagMaps.items():
        #     _d.update({_map.get(_pF, None): _map.get(_pT, None)})
        # # ensure tags not implemented will be mapped to None
        # if None in _d:
        #     _d.update({None:None})
        # # cls.print_cm_log("Mapping to xcTags", _d, level=3, depth=2)
        # if getAll:
        #     _d.pop(None, None)
        #     return tuple(_d.values())
        # if len(tags) == 1:
        #     return (_d.get(tags[0], None),)
        # return tuple(_d.get(_t, None) for _t in tags)
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