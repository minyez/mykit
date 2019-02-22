# -*- coding: utf-8 -*-
'''define classes and functions related to plane wave basis setup
'''

class pw:
    '''the class that defines variables that control the plane wave basis
    '''

    __encutPw = 0.0
    __unit = 'eV'

    def __init__(self, **kwargs):
        if "encutPw" in kwargs:
            self.__encutPw = kwargs["encutPW"]
        if "unit" in kwargs:
            self.__unitCutPw = kwargs["unit"]