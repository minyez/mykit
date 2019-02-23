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
        pass

    def parse_tags(self, **xctags):
        self.__parse_xctags(**xctags)

    def __parse_xctags(self, **xctags):
        pass

    @property
    def xcTags(self):
        return self.__xcTags