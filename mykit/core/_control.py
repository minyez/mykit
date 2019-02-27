# coding = utf8

from abc import ABC, abstractmethod
from mykit.core.log import verbose

class ControllerError(Exception):
    pass

class prog_mapper(ABC):
    '''Abstract base class for the tag and value mapping utility of controller classes.

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

    @abstractmethod
    def parse_tags(self, **kwargs):
        raise NotImplementedError

    @abstractmethod
    def pop_tags(self, *args):
        raise NotImplementedError

    @abstractmethod
    def delete_tags(self, *args):
        raise NotImplementedError

    @abstractmethod
    def _get_one_mykit_tag(self, tagName):
        raise NotImplementedError

    @abstractmethod
    def tag_vals(self, *args):
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


def check_valid_map(mapDict, basename=None):
    '''Check if a dict is a valid tag mapping object.

    Namely, a tag mapping object is a dict with mykit tag (str) and its mapping dict as key-value pair.
    When basename is specified, each mapping dict of tag 'm' must include a 'basename': 'm' member

    Args:
        mapDict (dict): a tag mapping object
        basename (str)

    Returns:
        True if mapDict is a valid mapping object, False otherwise
    '''
    try:
        assert isinstance(mapDict, dict)
        for k, v in mapDict.items():
            # key-value pair should be str:dict
            assert isinstance(k, str)
            assert isinstance(v, dict)
            if basename:
                assert basename in v.keys()
                assert k == v[basename]
    except AssertionError:
        return False
    else:
        return True


def tags_mapping(mapDict, progFrom, progTo, *tags, getAll=False):
    '''
    Args:
        m (dict): dictionary, with key-value pair of `mykit tag name: dict as tag map`
        progFrom (str)
        progTo (str)
    '''
    if len(tags) == 0 and not getAll:
        return tuple()
    if not check_valid_map(mapDict):
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


def parse_to_tagdict(tvDict, tagMaps, progName, **kwargs):
    '''

    Args:
        tvDict (dict) : the dict to parse tag-value pair into
        tagMaps (dict) : the mapping of mykit tags to those in program `progName`
        progName (str) : the program to which the keywords(tags) in kwargs should belong to
        kwargs : the tag-value to parse into
    '''
    try:
        assert isinstance(tvDict, dict)
    except AssertionError:
        raise ControllerError
    for _origTag, _v in kwargs.items():
        if _origTag == None:
            continue
        # check mykit tags
        elif _origTag in tagMaps.keys():
            tvDict.update({_origTag:_v})
        else:
        # check program-specific tags
        # ? Can be optimized
            for _t, _tmap in tagMaps.items():
                if _origTag == _tmap.get(progName, None):
                    tvDict.update({_t: _v})
                    break
    

def extract_from_tagdict(tvDict, tagMaps, progFrom, *tag, delete=False):
    pass


def vals_mapping(m, progFrom, progTo, *val, getAll=False):
    raise NotImplementedError