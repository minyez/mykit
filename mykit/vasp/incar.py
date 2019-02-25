# coding = utf-8
'''this module implements the ``incar`` class that manages the VASP input file INCAR
'''
import copy
import os
import re
from mykit.core.utils import check_duplicates_in_tag_tuple, trim_comment
from mykit.core.planewave import plane_wave_control, planewaveError
from mykit.core.xc import xc_control, xcError
# from mykit.core.program import program


class incarError(Exception):
    pass


# TODO check ion tags
class incar(plane_wave_control, xc_control):
    '''manage tags and IO of VASP input file INCAR
    '''
    # For VASP tags that is not easily to find analogy in other programs
    comment = ''
    vaspOnlyTags = ()
    tagElecBasic = (
                    "ISTART","LWAVE","IALGO","EDIFF","ENCUT","NBANDS",
                    "ICHARG","LCHARG","PREC","ALGO","ISPIN","NELM",
                    "ISMEAR","SIGMA","NELECT","LREAL",'NELMIN',
                    'MAGMOM','LMAXTAU',
                   )
    tagElecXc = (
                    "GGA","LHFCALC","PRECFOCK",'METAGGA',
                    'AEXX','AGGAX','AGGAC','ALDAC',
                )
    tagElecAdv = (
                    "ENCUTGW","LASPH","LPEAD","LOPTICS","ENCUTGWSOFT",
                    'NKRED','NKREDX','NKREDY','NKREDZ',"NOMEGA",
                    'NGX','NGY','NGZ','NGXF','NGYF','NGZF',
                    'ODDONLY','EVENONLY',"NBANDSGW","OMEGATL",
                    "OMEGAMAX","ANTIRES","NBANDSO","NBANDSV",
                 )
    tagIonBasic = (
                    'NSW','IBRION','ISIF','ISYM','EDIFFG','POTIM',
                  )
    tagIonSlab = (
                    "LDIPOL","DIPOL","IDIPOL",'LVHAR','LVTOT',
                 )
    tagPara = (
                    'NPAR','NCORE','KPAR',
              )
    tagNotCateg = (
                    'LMIXTAU','NWRITE','SYSTEM',
                    'AMIX','BMIX','IMIX','TIME','SMASS',
                  )
    __tagNotImple= ()
    tagAll = tagElecBasic + tagElecAdv + tagElecXc + \
        tagIonBasic + tagIonSlab + tagPara + tagNotCateg
    # ? Maybe move this duplicate check to unittest
    __hasDup = check_duplicates_in_tag_tuple(tagAll)
    if __hasDup > 0:
        raise incarError("Found tag duplicate in VASP tags. '{}' at Index {}".format(tagAll[__hasDup], __hasDup))
    # Map INCAR tags to tags of base classes
    tagMap2Base = {}
    __availMap2Pw = plane_wave_control.map2pwtags(*tagAll, progFrom="vasp")
    __availMap2Xc = xc_control.map2xctags(*tagAll, progFrom="vasp")
    for _m in [__availMap2Pw, __availMap2Xc]:
        for _i, _v in enumerate(_m):
            if _v != None:
                tagMap2Base.update({tagAll[_i]:_v})
    __vaspTags = {}

    # xc_dict = {
    #             # 'LDA'    : {'GGA'     : 'CA'} ,
    #             # 'PBE'    : {'GGA'     : 'PE'} ,
    #             # 'PBESOL' : {'GGA'     : 'PS'} ,
    #             # 'RPBE'   : {'GGA'     : 'RP'} ,
    #             # 'SCAN'   : {'METAGGA' : 'SCAN', 'LASPH' : '.TRUE.'} ,
    #             # 'PBE0'   : self.tag_PBE0  ,
    #             # 'HSE06'  : self.tag_HSE06 ,
    #             # 'HF'     : self.tag_HF
    #             }

    def __init__(self, **incarArgs):
        '''Initialize

        First, filter all incarargs to get all pwTags, xcTags, and remove their VASP equivalents,
        and then add those tags not implemented in the plane_wave tag mapping.
        '''
        super(incar, self).__init__("vasp", **incarArgs)
        self.parse_tags(**incarArgs)

    def __getattr__(self, attr):
        # pw and xc tags are not achievable by attribute
        _attr = attr.upper()
        if _attr not in self.tagAll:
            raise AttributeError("Invalid VASP tag: {}".format(_attr))
        return self.tag_vals(_attr)[0]
    
    def __setattr__(self, attr, value):
        # pw and xc tags are not achievable by attribute
        _attr = attr.upper()
        if _attr not in self.tagAll:
            self.print_warn("Invalid VASP tag: {}. Change nothing.".format(_attr),depth=0,level=1)
        else:
            self.parse_tags(**{_attr: value})
            

    def parse_tags(self, **keyval):
        if "comment" in keyval:
            self.comment = keyval.pop("comment")
        self.print_log("incar Parsing ", keyval, depth=0, level=3)
        plane_wave_control.parse_tags(self, "vasp", **keyval)
        xc_control.parse_tags(self, "vasp", **keyval)
        self.__parse_vasptags(**keyval)
        self.print_log("incar Finish parsing.\n", depth=0, level=3)
        # self.print_log("  pwTags",self.pwTags,depth=2,level=3)
        # self.print_log("  xcTags",self.xcTags,depth=2,level=3)
        # self.print_log("vaspTags",self.__vaspTags,depth=2,level=3)

    def __parse_vasptags(self, **tags):
        __tagFiltered = []
        if len(tags) != 0:
            __tagFiltered = self.__filter_tags_incar_not_pw_xc(*tags.keys())
        for _k in __tagFiltered:
            self.__vaspTags.update({_k:tags[_k]})
        self.print_log("End __parse_vasptags. vaspTags: ", self.__vaspTags, depth=1, level=3)

    def pop_tags(self, *tags):
        return self.__tag_vals(*tags, delete=True)

    def delete_tags(self, *tags):
        self.__tag_vals(*tags, delete=True)

    def tag_vals(self, *tags):
        return self.__tag_vals(*tags)

    # TODO implement delete
    def __tag_vals(self, *tags, delete=False):
        assert isinstance(delete, bool)
        self.print_log("In tag_vals (incar)", level=3, depth=0)
        if len(tags) == 0:
            return []
        __vals = [None,] * len(tags)
        if delete:
            __pwTagVals = plane_wave_control.pop_tags(self, "vasp", *tags)
            __xcTagVals = xc_control.pop_tags(self, "vasp", *tags)
        else:
            __pwTagVals = plane_wave_control.tag_vals(self, "vasp", *tags)
            __xcTagVals = xc_control.tag_vals(self, "vasp", *tags)
        __vaspTagVals = self.__vasptag_vals(*tags, delete=delete)
        __search = ( __pwTagVals, __xcTagVals, __vaspTagVals)
        self.print_log("Searching value in ", __search, depth=2, level=2)
        for _i, _v in enumerate(__vals):
            if _v is None:
                for _sList in __search:
                    if _sList[_i] != None:
                        __vals[_i] = _sList[_i]
                        break
        return __vals
                
    def __vasptag_vals(self, *tags, delete=False):
        _vals = []
        if len(tags) == 0:
            return _vals
        self.print_log("In __vasptag_vals, search {} tags of VASP: ".format(len(tags)), tags, level=3, depth=1)
        _vals = list(map(self.__vaspTags.get, tags))
        if delete:
            for _t in tags:
                self.__vaspTags.pop(_t, None)
        return _vals

    # @classmethod
    def __filter_tags_incar_not_pw_xc(self, *tags):
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
        __tagsDict = {}
        # Remove duplicate by constructing dict
        for _i,_t in enumerate(tags):
            __tagsDict.update({_t:_i})
        # Filter xc and pw tags that have INCAR map
        # INCAR tags that have common xc and pw tag mapping will be filtered as well
        # ? or filter all xc and pw tags?
        for _k, _v in __tagsDict.items():
            if _k in self.tagMap2Base.values() or _k in self.tagMap2Base.keys():
                pass
            elif _k in self.tagAll:
                _filter1.append(_k)
        # if len(_filter1) == 0:
        return _filter1
        # _filter2 = []
        # for _t in _filter1:
        #     if _t not in self.tagMap2Base:
        #         _filter2.append(_t)
        # return _filter2

    # TODO test write method
    def write(self, pathIncar="INCAR", backup=False, suffix="_bak"):
        '''Write to INCAR file at pathIncar. 

        Args:
            pathIncar (str) : the path to write the INCAR file.
            backup (bool) : when there is file 
            suffix (str) : suffix of the backup INCAR file
        '''
        _name = pathIncar
        try:
            assert not os.path.isdir(_name)
        except AssertionError:
            raise incarError("The path to write file is a directory.")
        if os.path.isfile(_name) and backup:
            _name = _name + suffix.strip()
        __all = self.tag_vals(*self.tagAll)
        with open(_name,'w') as f:
            if self.comment != '':
                print(self.comment, file=f)
            for _i, _v in enumerate(__all):
                if _v != None:
                    print(self.tagAll[_i], "=", _v, file=f)
            

    @classmethod
    def analyze_incar_line(cls, lineNum, incarLine, filePath=None):
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
        # 1. Trim line after # and !
        _l = trim_comment(_l, r"[\#\!]").strip()
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
                raise incarError("Bad tag pair on line {}. INCAR path: {}".format(lineNum+1, filePath))
            _k, _v = (_w.strip() for _w in kv)
            _k = _k.upper()
            # Check if tag is valid
            if _k not in cls.tagAll:
                raise incarError("Unknown INCAR tag on line {}: {}. INCAR path: {}".format(lineNum+1, _k, filePath))
            # Special case for SYSTEM key: the value is always string
            if _k == "SYSTEM":
                _kw.update({_k: _v})
                continue
            # Take care of the value. Raise for empty value
            if len(_v) == 0:
                raise incarError("Unknown INCAR tag on line {}. INCAR path: {}".format(lineNum+1, filePath))
            # Bool
            if _v.startswith('.'):
                if not _v.endswith('.'):
                    raise incarError("Bad bool tag on line {}. INCAR path: {}".format(lineNum+1, filePath))
                _v = _v[1:-1].lower()
                if _v == 'false':
                    _kw.update({_k:False})
                elif _v == 'true':
                    _kw.update({_k:True})
                else:
                    raise incarError("Bad bool tag on line {}. INCAR path: {}".format(lineNum+1, filePath))
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
                                    raise incarError("Bad multiplication tag on line {}. INCAR path: {}".format(lineNum+1, filePath))
                                try:
                                    _v = [float(_vList[1]),] * int(_vList[0])
                                except ValueError:
                                    raise incarError("Bad multiplication tag on line {}. INCAR path: {}".format(lineNum+1, filePath))
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
                                        raise incarError("Bad multiplication tag on line {}. INCAR path: {}".format(lineNum+1, filePath))
                                    try:
                                        _magSeg = [float(_vList[1]),] * int(_vList[0])
                                        newList.extend(_magSeg)
                                    except ValueError:
                                        raise incarError("Bad multiplication tag on line {}. INCAR path: {}".format(lineNum+1, filePath))
                    _kw.update({_k: newList})
        return _kw
                    
    @classmethod
    def read_from_file(cls, pathIncar="INCAR"):
        '''Generate ``incar`` instance from INCAR file
        '''
        __incarTags = {}
        if not os.path.isfile(pathIncar):
            raise incarError("Invalid path of INCAR file: {}".format(pathIncar))
        _f = open(pathIncar, 'r')
        _lines = _f.readlines()
        _f.close()
        for _i, _l in enumerate(_lines):
            __tags = cls.analyze_incar_line(_i, _l, filePath=pathIncar)
            for _k, _v in __tags.items():
                __incarTags.update({_k:_v})
        return cls(**__incarTags)
        
    @classmethod
    def minimal_electric(cls, pathIncar=None):
        pass

    @classmethod
    def minimal_ion_opt(cls, pathIncar=None):
        pass