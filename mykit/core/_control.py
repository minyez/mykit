# coding = utf8

from abc import ABC, abstractmethod
from mykit.core.log import verbose

class ControllerError(Exception):
    pass

class prog_mapper(ABC):
    '''Abstract base class for the mapping utility of controller classes.

    This ABC is aimed to provide the basic structure of any control class
    that is used to convert common input tags between programs.

    TODO:
        Values mapping
    '''
    @property
    @abstractmethod
    def _tagMaps(self):
        return NotImplementedError

    @classmethod
    @abstractmethod
    def map_tags(cls, *tags, progFrom="mykit", progTo="mykit", getAll=False):
        # return tags_mapping(cls._tagMaps, progFrom, progTo, *tags, getAll=False)
        raise NotImplementedError

    @classmethod
    def map_to_mykit_tags(cls, *tags, progFrom="mykit", getAll=False):
        '''Map some tags in a program to mykit tags, through map_tags method
        '''
        return cls.map_tags(*tags, progFrom=progFrom, progTo="mykit", getAll=getAll)
    
    @classmethod
    def map_from_mykit_tags(cls, *tags, progTo="mykit", getAll=False):
        '''Map some program-specific tags to mykit tags, through map_tags method
        '''
        return cls.map_tags(*tags, progFrom="mykit", progTo=progTo, getAll=getAll)

# TODO parser ABC for base classes
class parser(ABC):
    '''Abstract base class for the value parser utility of controller classes.
    '''
    @abstractmethod
    def parse_tags(self, *kwargs):
        raise NotImplementedError

    @abstractmethod
    def pop_tags(self, *kwargs):
        raise NotImplementedError

    @abstractmethod
    def delete_tags(self, *kwargs):
        raise NotImplementedError

    @abstractmethod
    def tag_vals(self, *kwargs):
        raise NotImplementedError


def tags_mapping(mapDict, progFrom, progTo, *tags, getAll=False):
    '''
    Args:
        m (dict): dictionary, with key-value pair of `mykit tag name: dict as tag map`
        progFrom (str)
        progTo (str)
    '''
    if len(tags) == 0 and not getAll:
        return tuple()
    try:
        assert isinstance(mapDict, dict)
        for _v in mapDict.values():
            assert isinstance(_v, dict)
    except AssertionError:
        raise ControllerError
    verbose.print_cm_log("Use mapDict:", mapDict, level=4, depth=1)
    verbose.print_cm_log("To map tags:", tags, level=4, depth=1)
    _pF = progFrom
    _pT = progTo
    _d = {}
    for _map in mapDict.values():
        _d.update({_map.get(_pF, None): _map.get(_pT, None)})
    if None in _d:
        _d.update({None:None})
    if getAll:
        _d.pop(None, None)
        return tuple(_d.values())
    if len(tags) == 1:
        return (_d.get(tags[0], None),)
    return tuple(_d.get(_t, None) for _t in tags)


def vals_mapping(m, progFrom, progTo, *val, getAll=False):
    raise NotImplementedError