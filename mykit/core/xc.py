# coding = utf8
'''
'''

class xcError(Exception):
    pass


class xc_control:
    '''base class that controls exchange-correlation setup
    '''

    __xcTagMaps = {
                    "gga" : {"n a":"gga", "vasp":"GGA"},
                    "metagga" : {"n a":"metagga", "vasp":"METAGGA"}
                  }
    __xcTags = {}

    def __init__(self, **xctags):
        self.__parse_xctags(**xctags)

    def parse_tags(self, **xctags):
        self.__parse_xctags(**xctags)

    def __parse_xctags(self, **xctags):
        pass

    def tag_vals(self, *tags, progName="n a"):
        raise NotImplementedError
        # return self.__xctag_vals(*tags, progName=progName)

    def __xctag_vals(self, *tags, progName="n a"):
        raise NotImplementedError

    @classmethod
    def map_tags_in_xc(cls, *tags, progFrom="n a", progTo="n a", getAll=False):
        raise NotImplementedError
    
    @classmethod
    def map2xctags(cls, *tags, progFrom="n a", getAll=False):
        cls.map_tags_in_xc(cls, *tags, progFrom=progFrom, progTo="n a", getAll=getAll)

    @classmethod
    def map_from_xctags(cls, *tags, progTo="n a", getAll=False):
        cls.map_tags_in_xc(cls, *tags, progTo=progTo, progFrom="n a", getAll=getAll)

    @property
    def xcTags(self):
        return self.__xcTags