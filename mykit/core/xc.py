# coding = utf8
'''
'''

class xcError(Exception):
    pass


class xc_control:
    '''base class that controls exchange-correlation setup
    '''

    __xcTagMaps = {}
    __xcTags = {}

    def __init__(self, **xctags):
        pass