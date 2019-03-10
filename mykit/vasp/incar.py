# coding = utf-8
'''this module implements the ``incar`` class that manages the VASP input file INCAR
'''
import copy
import json
import os
import re

from mykit.core.ion import ion_control
from mykit.core.planewave import planewave_control
from mykit.core.utils import check_duplicates_in_tag_tuple, trim_comment
from mykit.core.xc import xc_control

# from mykit.core.program import program

_incar_controllers = (
    planewave_control, 
    xc_control, 
    ion_control,
    )

class IncarError(Exception):
    pass


class incar(*_incar_controllers):
    '''manage tags and IO of VASP input file INCAR
    '''
    # For VASP tags that is not easily to find analogy in other programs
    _metadata = os.path.join(os.path.dirname(__file__), "metadata", "incartags.json")
    with open(_metadata, 'r') as h:
        _tagdoc = json.load(h)
    tagAll = []
    docAll = []
    for _t in _tagdoc.values():
        tagAll.extend(_t.keys())
        docAll.extend(_t.values())
    # ? Maybe move this duplicate check to unittest
    _hasDup = check_duplicates_in_tag_tuple(tagAll)
    if _hasDup > 0:
        raise IncarError("Found tag duplicate in VASP tags. '{}' at Index {}".format(tagAll[_hasDup], _hasDup))
    # Map INCAR tags to tags of base classes
    incar2mykit = {}
    for _c in _incar_controllers:
        _m = _c.map_to_mykit_tags(*tagAll, progFrom="vasp")
        for _i, _v in enumerate(_m):
            if _v != None:
                incar2mykit.update({tagAll[_i]:_v})

    def __init__(self, **incarArgs):
        '''Initialize

        First, filter all incarargs to get all pwTags, xcTags, and remove their VASP equivalents,
        and then add those tags not implemented in the controller tag mapping.
        '''
        self.comment = ''
        self._vaspTags = {}
        for _c in _incar_controllers:
            _c.__init__(self, "vasp", **incarArgs)
        self.parse_tags(**incarArgs)

    @property
    def tags(self):
        '''get all parsed tags
        '''
        _ret = []
        __alltv = self.tag_vals(*self.tagAll)
        for tv in zip(self.tagAll, __alltv):
            if tv[1] != None:
                _ret.append(tv)
        return dict(_ret)

    def __str__(self):
        return "{}".format(self.tags)

    def __getitem__(self, tag):
        try:
            assert tag in self.tagAll
        except AssertionError:
            raise IncarError("Invalid tag for INCAR to extract: {}".format(tag))
        return self._tag_vals(tag)[0]

    def __setitem__(self, tag, val):
        try:
            assert tag in self.tagAll
        except AssertionError:
            raise IncarError("Invalid tag for INCAR to parse: {}".format(tag))
        self.parse_tags(**{tag: val})

    def __repr__(self):
        return self.__str__()

    def parse_tags(self, **keyval):
        # print(keyval)
        if "comment" in keyval:
            self.comment = keyval.pop("comment")
        self.print_log("incar Parsing ", keyval, depth=0, level=3)
        for _c in _incar_controllers:
            _c.parse_tags(self, "vasp", **keyval)
        self._parse_vasptags(**keyval)
        self.print_log("incar Finish parsing.\n", depth=0, level=3)
        # self.print_log("  pwTags",self.pwTags,depth=2,level=3)
        # self.print_log("  xcTags",self.xcTags,depth=2,level=3)
        # self.print_log("vaspTags",self._vaspTags,depth=2,level=3)

    def _parse_vasptags(self, **tags):
        __tagFiltered = []
        if len(tags) != 0:
            __tagFiltered = self._filter_tags_incar_not_mykit(*tags.keys())
        for _k in __tagFiltered:
            self._vaspTags.update({_k:tags[_k]})
        # self.print_log("End _parse_vasptags. vaspTags: ", self._vaspTags, depth=1, level=3)

    def pop_tags(self, *tags):
        return self._tag_vals(*tags, delete=True)

    def delete_tags(self, *tags):
        self._tag_vals(*tags, delete=True)

    def tag_vals(self, *tags):
        return self._tag_vals(*tags)

    def _tag_vals(self, *tags, delete=False):
        assert isinstance(delete, bool)
        self.print_log("In tag_vals (incar)", level=3, depth=0)
        if len(tags) == 0:
            return []
        vals = [None,] * len(tags)
        _search = []
        for _c in _incar_controllers:
            if delete:
                _search.append(_c.pop_tags(self, "vasp", *tags))
            else:
                _search.append(_c.tag_vals(self, "vasp", *tags))
        _vaspTagVals = self._vasptag_vals(*tags, delete=delete)
        _search.append(_vaspTagVals)
        # self.print_log("Searching value in ", _search, depth=2, level=2)
        for _i, _v in enumerate(vals):
            if _v is None:
                for _sList in _search:
                    if _sList[_i] != None:
                        vals[_i] = _sList[_i]
                        break
        return vals
                
    def _vasptag_vals(self, *tags, delete=False):
        _vals = []
        if len(tags) == 0:
            return _vals
        # self.print_log("In _vasptag_vals, search {} tags of VASP: ".format(len(tags)), tags, level=3, depth=1)
        _vals = list(map(self._vaspTags.get, tags))
        if delete:
            for _t in tags:
                self._vaspTags.pop(_t, None)
        return _vals

    def _filter_tags_incar_not_mykit(self, *tags):
        '''Get VASP tags that are implemented in incar and have no correspondents in base class tags.
        
        Common plane_wave and xc tags, and their VASP correspondent,
        will be filtered out.

        Returns:
            list
        
        Note:
            Order is not kept.
        '''
        _filter1 = []
        if len(tags) == 0:
            return _filter1
        tagsDict = {}
        # Remove duplicate by constructing dict
        for _i,_t in enumerate(tags):
            tagsDict.update({_t:_i})
        # Filter xc and pw tags that have INCAR map
        # INCAR tags that have common xc and pw tag mapping will be filtered as well
        # ? or filter all xc and pw tags?
        for _k, _v in tagsDict.items():
            if _k in self.incar2mykit.values() or _k in self.incar2mykit.keys():
                pass
            elif _k in self.tagAll:
                _filter1.append(_k)
        return _filter1

    def __print(self, f):
        '''Print the incar tag-value dict to file-type f in the INCAR format

        Args:
            f (file-type)
        '''
        # __all = self.tag_vals(*self.tagAll)
        if self.comment != '':
            print(self.comment, file=f)
        # for _i, _v in enumerate(__all):
        for k, v in self.tags.items():
            if v != None:
                if isinstance(v, (list, tuple)):
                    print(k, "=", *v, file=f)
                elif isinstance(v, dict):
                    print(k, "=", **v, file=f)
                elif isinstance(v, bool):
                    print(k, "=", {True:".TRUE.", False:".FALSE."}[v], file=f)
                else:
                    print(k, "=", v, file=f)

    def print(self):
        '''Preview the INCAR output
        '''
        from sys import stdout
        self.__print(stdout)

    # TODO test write method, i.e. __print method
    def write(self, pathIncar="INCAR", backup=False, suffix="_bak"):
        '''Write to INCAR file at pathIncar. 

        Args:
            pathIncar (str) : the path to write the INCAR file.
            backup (bool) : when set True and file ``pathIncar`` exists, the original file will
                be backuped as ``pathIncar``+``suffix``
            suffix (str) : suffix of the backup INCAR file
        '''
        _name = pathIncar
        try:
            assert not os.path.isdir(_name)
        except AssertionError:
            raise IncarError("The path to write INCAR is a directory.")
        if os.path.isfile(_name) and backup:
            _bakname = _name + suffix.strip()
            os.rename(_name, _bakname)
        with open(_name, 'w') as f:
            self.__print(f)
            
    @classmethod
    def analyze_incar_line(cls, incarLine, lineNum=-1, filePath=None):
        '''Analyze one INCAR line
        '''
        _kw = {}
        _l = incarLine.strip()
        if len(_l) == 0:
            return {}
        else:
            if lineNum == 0:
                # If '=' is not found in the first line, set it as the comment.
                if '=' not in _l:
                    return {"comment":_l}
                # Otherwise, check if the first word before "=" belongs to a known VASP tag in tagAll
                # if not, treat the line as a comment as well
                else:
                    _w = _l.split('=')[0]
                    if _w not in cls.tagAll:
                        return {"comment":_l}
            else:
                # if '=' is not found in lines after first, skip it
                if '=' not in _l:
                    return {}
        # Here start treat a line that possibly contain INCAR tags
        # 1. Trim line after #, ! and (
        _l = trim_comment(_l, r"[\#\!\(]").strip()
        if len(_l) == 0:
            return _kw
        # 2. Split the line by semicolon.
        _words = _l.split(';')
        # 3. iterate for each segment
        for _w in _words:
            if _w == '':
                continue
            kv = _w.split("=")
            # At this stage, this must be 2-member pair as there can be one equality only.
            if len(kv) != 2:
                raise IncarError("Bad tag pair on line {}. INCAR path: {}".format(lineNum+1, filePath))
            _k, _v = (_w.strip() for _w in kv)
            _k = _k.upper()
            # Check if tag is valid
            if _k not in cls.tagAll:
                raise IncarError("Unknown INCAR tag on line {}: {}. INCAR path: {}".format(lineNum+1, _k, filePath))
            # Special case for SYSTEM key: the value is always string
            if _k == "SYSTEM":
                _kw.update({_k: _v})
                continue
            # Take care of the value. Raise for empty value
            if len(_v) == 0:
                raise IncarError("Unknown INCAR tag on line {}. INCAR path: {}".format(lineNum+1, filePath))
            # Bool
            if _v.startswith('.'):
                if not _v.endswith('.'):
                    raise IncarError("Bad bool tag on line {}. INCAR path: {}".format(lineNum+1, filePath))
                _v = _v[1:-1].lower()
                if _v == 'false':
                    _kw.update({_k:False})
                elif _v == 'true':
                    _kw.update({_k:True})
                else:
                    raise IncarError("Bad bool tag on line {}. INCAR path: {}".format(lineNum+1, filePath))
            else:
                # TODO this part can be extracted to form an independent string function
                # check if _v is a list of value
                _vList = _v.split()
                if len(_vList) == 1:
                    try:
                        # finally a integer
                        _v = int(_v)
                    except ValueError:
                        try:
                            # finally a float
                            _v = float(_v)
                        except ValueError:
                            _vList = _v.split("*")
                            if len(_vList) == 1:
                                # finally a string
                                pass
                            else:
                                # a single string with * in it, such as MAGMOM
                                if len(_vList) != 2:
                                    raise IncarError("Bad multiplication tag on line {}. INCAR path: {}".format(lineNum+1, filePath))
                                try:
                                    _v = [float(_vList[1]),] * int(_vList[0])
                                except ValueError:
                                    raise IncarError("Bad multiplication tag on line {}. INCAR path: {}".format(lineNum+1, filePath))
                    _kw.update({_k:_v})
                else:
                    newList = []
                    # _break = False
                    for _seg in _vList:
                        if _seg == ',':
                            # _break = True
                            continue
                        try:
                            # finally a integer
                            _intSeg = int(_seg)
                            newList.append(_intSeg)
                        except ValueError:
                            try:
                                # finally a float
                                _floatSeg = float(_seg)
                                newList.append(_floatSeg)
                            except ValueError:
                                _segSplit = _seg.split("*")
                                if len(_segSplit) == 1:
                                    # finally a string
                                    newList.append(_seg)
                                else:
                                    # a string with * in it, such as MAGMOM: 48*5.0 16*1.0
                                    if len(_vList) != 2:
                                        raise IncarError("Bad multiplication tag on line {}. INCAR path: {}".format(lineNum+1, filePath))
                                    try:
                                        _magSeg = [float(_vList[1]),] * int(_vList[0])
                                        newList.extend(_magSeg)
                                    except ValueError:
                                        raise IncarError("Bad multiplication tag on line {}. INCAR path: {}".format(lineNum+1, filePath))
                    _kw.update({_k: newList})
        return _kw
                    
    # * Factory methods
    @classmethod
    def read_from_file(cls, pathIncar="INCAR"):
        '''Generate ``incar`` instance from INCAR file
        '''
        __incarTags = {}
        if not os.path.isfile(pathIncar):
            raise IncarError("Invalid path of INCAR file: {}".format(pathIncar))
        _f = open(pathIncar, 'r')
        _lines = _f.readlines()
        _f.close()
        for _i, _l in enumerate(_lines):
            __tags = cls.analyze_incar_line(_l, lineNum=_i, filePath=pathIncar)
            for _k, _v in __tags.items():
                __incarTags.update({_k:_v})
        return cls(**__incarTags)
        
    @classmethod
    def minimal_scf(cls, xc=None, nproc=1, **kwargs):
        '''Create an ``incar`` instance with minimal reasonable tags for SCF

        Args:
            xc (str)
            nproc (int): number of processors to use
            kwargs: other tags you want to set manually.
                 It is not recommended as it will overwrite the default tags, and
                 in this case you may consider other factory methods.

        Returns:
            incar object
        '''
        _scftags = {
            "LREAL": False,
            "EDIFF": 1e-6,
            "PREC": "Normal",
            "ISTART": 0,
            "ICHARG": 2,
            "LMAXMIX": 6,
            }
        _scftags.update(kwargs)
        if xc != None:
            _xctags = _get_xc_tags_from_xcname(xc)
            _scftags.update(_xctags)
        _paratags = _get_para_tags_from_nproc(nproc)
        _scftags.update(_paratags)
        return cls(**_scftags)

    @classmethod
    def minimal_ion_opt(cls, **kwargs):
        pass


def _get_xc_tags_from_xcname(xc, ignore_error=True):
    '''Get necessary VASP tags for running calculation of ``xc`` functional

    Currently support xc name

        local: LDA, PBE, PBEsol, RPBE
        meta-GGA: SCAN
        hybrid: PBE0, HSE06, HF

    Args:
        xc (str) : the name of exchange-correlation functional, case-insensitive.
        ignore_error (bool): if raise ``IncarError`` when the xc name is not supported.
            If set False and ``xc`` is not supported, an empty dict will be returned.
        
    Returns:
        dict
    '''
    _tag_PBE0 = (("GGA", "PE"), ("ALGO", "ALL"), ("LHFCALC", True), ("TIME", 0.4), ("PRECFOCK", "Normal"))
    _tags = {
        "LDA": {"GGA": "CA"}, 
        "PBE": {"GGA": "PE"},
        "PBESOL": {"GGA": "PE"},
        "RPBE": {"GGA": "RP"},
        "SCAN": {"METAGGA": "SCAN", "LASPH": True}, 
        "PBE0": dict(_tag_PBE0),
        "HSE06": dict(_tag_PBE0 + (("HFSCREEN", 0.2),)),
        "HF": dict(_tag_PBE0 + (("AEXX", 1.0),)),
        }
    _xc = xc.upper()
    if ignore_error:
        return _tags.get(_xc, {})
    else:
        raise IncarError("XC name is not supported: {}".format(xc))


def _get_para_tags_from_nproc(nproc):
    '''Get the value of parallelization related tags from number of processors.

    The returned dict will include two keys, KPAR and NPAR, if ``nproc`` is not prime.
    Their product equals ``nproc``. NPAR will be no smaller than KPAR.
    Otherwise, i.e. ``nproc`` is prime, an empty dict will be returned.

    Args:
        nproc (int): number of processors to use

    Returns:
        dict
    '''
    from math import sqrt
    assert isinstance(nproc, int)
    assert nproc > 0
    _paratags = {}
    if nproc > 1:
        _kpar = int(sqrt(nproc))
        while _kpar > 1:
            if nproc%_kpar == 0:
                _paratags.update({"KPAR": _kpar, "NPAR": int(nproc/_kpar)})
                break
            _kpar -= 1
    return _paratags
