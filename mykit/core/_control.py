# coding = utf8

from abc import ABC, abstractmethod

from mykit.core.log import Verbose


class ControllerError(Exception):
    pass

class prog_mapper(ABC):
    '''Abstract base class for the tag and value mapping utility of controller classes.

    This ABC is aimed to provide the basic structure of any control class
    that is used to convert common input tags between programs.

    Notes:
        A `_tagMaps` abstract property is defined and should be overwriten by tag mapping
        object when defining a new controller class. 
        However, the class should not directly use it to parse and extract tags,
        as some program input may need to inherit from several controller classes,
        then when inheriting, duplicate _tagMaps will lose its reference and cause failure
        in mapping of tags in some controller.

        It is solved by referring to `_tagMaps` by some other internal variable,
        e.g. `_xcTagMaps` in `xc_control`, and use it when calling, e.g. mapping functions.

        This issue also appears for methods and property, but not classmethod.

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
    def parse_tags(self, progName, **kwargs):
        raise NotImplementedError

    @abstractmethod
    def pop_tags(self, progName, *args):
        raise NotImplementedError

    @abstractmethod
    def delete_tags(self, progName, *args):
        raise NotImplementedError

    @abstractmethod
    def _get_one_mykit_tag(self, tagName):
        raise NotImplementedError

    @abstractmethod
    def tag_vals(self, progName, *args):
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

    Namely, a tag mapping object is a dict with string keyword and a dict as key-value pair.
    When basename is specified, value (dict) of keyword 'm' must include {'basename': 'm'}

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
    Verbose.print_cm_log("Use mapDict:", mapDict, level=4, depth=1)
    Verbose.print_cm_log("To map tags:", tags, level=4, depth=1)
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


def parse_to_tagdict(tvDict, tagMaps, progName, **tagvals):
    '''Parse tag values to the dict, according to tagMaps.

    Note:
        If the mykit tag and its program-specific equivalent are parsed at
        the same time, the program-specific tag value will be preferred.

    Args:
        tvDict (dict) : the dict to parse tag-value pair into
        tagMaps (dict) : the tag mapping object
        progName (str) : the program to which the keywords(tags) in kwargs should belong to
        kwargs : the tag-value to parse into
    '''
    _parsed_mykit_tag = []
    try:
        assert isinstance(tvDict, dict)
    except AssertionError:
        raise ControllerError
    for _origTag, _v in tagvals.items():
        if _origTag == None:
            continue
        # check if mykit tags. Parsed it only when its program-specific tag is not parsed at the same time
        elif _origTag in tagMaps.keys():
            if _origTag not in _parsed_mykit_tag:
                tvDict.update({_origTag:_v})
        else:
        # check program-specific tags
        # ? Can be optimized
            for _t, _tmap in tagMaps.items():
                if _origTag == _tmap.get(progName, None):
                    _parsed_mykit_tag.append(_t)
                    tvDict.update({_t: _v})
                    break
    

def extract_from_tagdict(controlClass, tvDict, progName, *tags, delete=False):
    '''Extract values from dict that stores mykit tag value by tag names of a program

    `controlClass` is required, as its class method map_to_mykit_tags will be used
    It also extract the value of mykit tag name, not only tags of progName.

    Args:
        controlClass (class) : the name of class which is a subclass of prog_mapper.
        tvDict (dict): the dict to manipulate which contains mykit tag name and value pair
        progName (str): the name of program to which the tags argument belong.
        tags (str): the names of tags to request
        delete (bool): if set True, the mykit tag will be deleted from tvDict
    '''
    try:
        assert issubclass(controlClass, prog_mapper)
    except:
        raise ControllerError("{}: Not a controller class, required as subclass of prog_mapper.".format(controlClass))
    if len(tags) == 0:
        return []
    mytags = controlClass.map_to_mykit_tags(*tags, progFrom=progName)
    # Verbose.print_cm_log("extract mykit tags: ", _xctags, level=3, depth=2)
    myvals = list(map(tvDict.get, mytags))
    # Verbose.print_cm_log("      their values: ", _vals, level=3, depth=2)
    for i, v in enumerate(myvals):
        if v == None:
            if tags[i] in tvDict:
                myvals[i] = tvDict[tags[i]]
    if delete:
        for i, v in enumerate(mytags):
            if v in tvDict:
                del tvDict[v]
            # in case there is mykit tag name in tags, delete it as well
            if tags[i] in tvDict:
                del tvDict[tags[i]]
    return myvals


def vals_mapping(m, progFrom, progTo, *val, getAll=False):
    raise NotImplementedError


def build_tag_map_obj(metadataFile, basename, filetype):
    '''Build a tag mapping object from a metadata file
    '''
    if filetype == "json":
        return _build_tag_map_obj_by_json(metadataFile, basename)
    else:
        raise NotImplementedError


def _build_tag_map_obj_by_json(jsonFile, basename):
    '''Build a tag mapping object from a JSON file

    Structure of data file:
        tag of basename: {"vasp": vasptag, "qe": qetag}
    '''
    import json

    _tagMaps = {}
    try:
        with open(jsonFile, 'r') as f:
            j = json.load(f)
    except json.JSONDecodeError:
        raise ControllerError("Fail to generate tag mapping from JSON file: {}".format(jsonFile))
    
    for k, v in j.items():
        try:
            assert isinstance(v, dict)
        except AssertionError:
            raise ControllerError("Bad dict in key {} of JSON file: {}".format(k, jsonFile))
        v.update({basename: k})
        _tagMaps.update({k: v})
    check_valid_map(_tagMaps, basename=basename)
    return _tagMaps
