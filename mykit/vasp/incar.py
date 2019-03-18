# coding = utf-8
'''this module implements the ``incar`` class that manages the VASP input file INCAR
'''
import copy
import json
import os
import re

from mykit.core.ion import ion_control
from mykit.core.planewave import planewave_control
from mykit.core.utils import check_duplicates_in_tag_tuple, trim_after
from mykit.core.xc import xc_control

# from mykit.core.program import program

_incar_controllers = (
    planewave_control, 
    xc_control, 
    ion_control,
    )


class IncarError(Exception):
    pass


class Incar(*_incar_controllers):
    '''manage tags and IO of VASP input file INCAR
    '''
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
        and then add to _vaspTags those tags not implemented in the controller tag mapping.
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

    def __repr__(self):
        return "{}".format(self.tags)

    def __str__(self):
        # __all = self.tag_vals(*self.tagAll)
        _ret = []
        if self.comment != '':
            _ret.append(self.comment)
        for k, v in self.tags.items():
            if v != None:
                _str = k + " = "
                if isinstance(v, (list, tuple)):
                    _str += ' '.join([str(x) for x in v])
                # elif isinstance(v, dict):
                #     print(k, "=", **v, file=f)
                elif isinstance(v, bool):
                    _str += {True:".TRUE.", False:".FALSE."}[v]
                else:
                    _str += str(v)
                _ret.append(_str)
        return '\n'.join(_ret)

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
        '''Get VASP tags that are implemented in incar but have no correspondents in base class tags.
        
        Common planewave, ion and xc tags, and their VASP correspondent,
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
            print(self.__str__(), file=f)
            
    @classmethod
    def analyze_incar_line(cls, incarLine, lineNum=-1, filePath=None):
        '''Analyze one INCAR line

        The steps follow:
            1. Trim line after #, ! and (
            2. Split the line by semicolon.
            3. iterate for each segment, each giving a tag-value pair
        '''
        _kw = {}
        _l = incarLine.strip()
        if len(_l) == 0:
            return {}
        else:
            if lineNum == 0:
                if '=' not in _l:
                    return {"comment":_l}
                else:
                    _w = _l.split('=')[0].strip()
                    if _w not in cls.tagAll:
                        return {"comment":_l}
            else:
                if '=' not in _l:
                    return {}
        # Here start treat a line that possibly contain INCAR tags
        _l = trim_after(_l, r"[\#\!\(]").strip()
        if len(_l) == 0:
            return _kw
        _words = _l.split(';')
        for _w in _words:
            if _w == '':
                continue
            kv = _w.split("=")
            # At this stage, this must be 2-member pair as there can be one equality only.
            if len(kv) != 2:
                raise IncarError("Bad tag pair on line {}. INCAR path: {}".format(lineNum+1, filePath))
            _k, _v = (_w.strip() for _w in kv)
            _k = _k.upper()
            if _k not in cls.tagAll:
                raise IncarError("Unknown INCAR tag on line {}: {}. INCAR path: {}".format(lineNum+1, _k, filePath))
            # Special case for SYSTEM key: the value is always string
            if _k == "SYSTEM":
                _kw.update({_k: _v})
                continue
            decoded = decode_incar_value(_v)
            if decoded is None:
                raise IncarError("Bad tag value on line {}. INCAR path: {}".format(lineNum+1, filePath))
            _kw.update({_k: decoded})
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
            __incarTags.update(__tags)
        return cls(**__incarTags)
        
    @classmethod
    def minimal_incar(cls, task, xc=None, nproc=1, **kwargs):
        '''Create an ``incar`` instance with minimal reasonable tags for particular task

        Args:
            task (str): the type of task, "band", "dos", "scf", "opt"
            xc (str): the type of xc.
            nproc (int): number of processors to use
            kwargs: other tags you want to set manually.
                Note that this will overwrite built-in defaults.

        Returns:
            incar object
        '''
        _taskDict = {
                "scf": (_get_minimal_scf_tags, ),
                "opt": (_get_minimal_scf_tags, _get_minimal_ion_opt_tags, ),
                "band": (_get_minimal_scf_tags, _get_minimal_band_tags, ),
                "dos": (_get_minimal_scf_tags, _get_minimal_dos_tags, ),
                "diag": (_get_minimal_scf_tags, _get_minimal_diag_tags, ),
            }
        minimals = _taskDict.get(task, None)
        if minimals is None:
            raise IncarError("Task name not supported: {}. Should be in {}".format(task, _taskDict.keys()))

        _tags = {}
        for minimal in minimals:
            _tags.update(minimal())
        if xc != None:
            _xctags = _get_xc_tags_from_xcname(xc)
            _tags.update(_xctags)
        _paratags = _get_para_tags_from_nproc(nproc)
        _tags.update(_paratags)
        _tags.update(kwargs)
        return cls(**_tags)


def _get_minimal_scf_tags():
    _scftags = {
        "LREAL": False,
        "EDIFF": 1e-6,
        "PREC": "Normal",
        "LMAXMIX": 6,
        }
    return _scftags


def _get_minimal_diag_tags():
    _diatags = {
        "NELM": 1,
        "ALGO": "Exact",
    }
    return _diatags


def _get_minimal_dos_tags():
    _dostags = {
        "ICHARG": 11,
        "LORBIT": 11,
        "ISTART": None,
        "NEDOS": 3000,
    }
    return _dostags


def _get_minimal_band_tags():
    _bandtags = {
        "ICHARG": 11,
        "LORBIT": 11,
        "ISTART": None,
    }
    return _bandtags


def _get_minimal_ion_opt_tags():
    _iontags = {
        "IBRION": 1,
        "EDIFFG": -0.01,
        "ISIF": 3,
        "NSW": 40,
        "MAXMIX": 80,
    }
    return _iontags


def _get_minimal_gw_tags():
    raise NotImplementedError


def _get_xc_tags_from_xcname(xc=None, ignore_error=True):
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
        "PBESOL": {"GGA": "PS"},
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


def decode_incar_value(vstr):
    '''decode a string of INCAR tag value into right data type

    Args:
        vstr (str): value string to decode
    '''
    boolDict = {'false': False, 'true': True}
    def __convert(v):
        # Bool
        if v.startswith('.'):
            if v.endswith('.'):
                b = vs[1:-1].lower()
                return boolDict.get(b, None)
        try:
            return int(v)
        except ValueError:
            try:
                return float(v)
            except ValueError:
                _vList = v.split("*")
                # String
                if len(_vList) == 1:
                    return v
                elif len(_vList) == 2:
                    # a single string with * in it, such as MAGMOM
                    try:
                        return [float(_vList[1]),] * int(_vList[0])
                    except ValueError:
                        pass
        return None

    vs = vstr.strip()
    value = None
    if len(vs) == 0:
        pass
    else:
        # check if string consists a list of value
        _vList = vs.split()
        if len(_vList) == 1:
            _v = _vList[0]
            value = __convert(_v)
        else:
            newList = []
            for seg in _vList:
                if seg == ',':
                    continue
                vseg = __convert(seg)
                if vseg is None:
                    return None
                elif isinstance(vseg, list):
                    newList.extend(vseg)
                else:
                    newList.append(vseg)
            value = newList
    return value
