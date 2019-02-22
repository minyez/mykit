# -*- coding: utf-8 -*-
'''Defined functions for printing mykit information
'''
from mykit.core.config import global_config
import sys

class verbose:
    '''Class that controls the level of information print.

    Warning and information have their own level of verbose.
    The levels can be set in the configure file, ``verbLog`` and ``verbWarn``,
    and also be adjusted by calling ``warnLevel`` and ``logLevel`` property
    
    level 0: Minimal
    level 1: Normal
    level 2: Check work
    level 3: Debug code

    Level set higher than the above number in configure file is not necessary.
    Level set higher than the above number by calling property will be ignored.
    '''

    _verbWarn, _verbLog, _indent = global_config.get('verbWarn', 'verbLog', 'logIndent')
    _indent = ' ' * _indent
    _prefix = {"log": '', "warn": ' WARNING!!!'}
    # range of verbose level
    __rangeW = range(0, 4)
    __rangeL = range(0, 4)

    def __init__(self):
        pass

    @property
    def warnLevel(self):
        return self._verbWarn
    @warnLevel.setter
    def warnLevel(self, value):
        if value in self.__rangeW:
            self._verbWarn = value

    @property
    def logLevel(self):
        return self._verbLog
    @logLevel.setter
    def logLevel(self, value):
        if value in self.__rangeL:
            self._verbLog = value

    def __print(self, *strs, printType="log", indentLevel=0, verbLevel=0, file=sys.stdout):
        assert isinstance(indentLevel, int)
        _pt = printType.lower()
        assert _pt in ["log", "warn"]
        _verbDict = {"log": self._verbLog, "warn": self._verbWarn}
        if verbLevel <= _verbDict[_pt]:
            _strPref = self._indent * indentLevel + self._prefix[_pt]
            print(_strPref, *strs, file=file)

    def print_log(self, *logStr, depth=0, level=0, file=sys.stdout):
        '''Print out log information with ``level`` larger than verbose setup ``logLevel``.
        '''
        self.__print(*logStr, printType="log", \
            indentLevel=depth, verbLevel=level, file=file)

    def print_warn(self, *warnStr, depth=0, level=0, file=sys.stdout):
        '''print out warning information with ``level`` larger than verbose setup ``warnLevel``.
        '''
        self.__print(*warnStr, printType="warn", \
            indentLevel=depth, verbLevel=level, file=file)
    
    @classmethod
    def print_cm_log(cls, *strs, indentLevel=0, verbLevel=0, file=sys.stdout):
        '''print log information for classmethod

        This helps the classmethod of the subclasses of verbose to print hierachal log.
        '''
        assert isinstance(indentLevel, int)
        _strPref = cls._indent * indentLevel + cls._prefix["log"]
        if verbLevel <= cls._verbLog:
            print(_strPref, *strs, file=file)

    @classmethod
    def print_cm_warn(cls, *strs, indentLevel=0, verbLevel=0, file=sys.stdout):
        '''print warning information for classmethod

        This helps the classmethod of the subclasses of verbose to print hierachal log.
        '''
        assert isinstance(indentLevel, int)
        _strPref = cls._indent * indentLevel + cls._prefix["warn"]
        if verbLevel <= cls._verbWarn:
            print(_strPref, *strs, file=file)
     