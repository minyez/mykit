# coding = utf-8
'''this module implements the ``incar`` class that manages the VASP input file INCAR
'''
import copy
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
    vaspOnlyTags = ()
    tagElecBasic = (
                    "ISTART","LWAVE","IALGO","EDIFF","ENCUT","NBANDS",
                    "ICHARG","LCHARG","PREC","ALGO","ISPIN","NELM",
                    "ISMEAR","SIGMA","NELECT","ENCUT","LREAL",'NELMIN',
                    'MAGMOM','LMAXTAU',
                   )
    tagElecXc = (
                    "GGA","LHFCALC","PRECFOCK",'METAGGA',
                    'AEXX','AGGAX','AGGAC','ALDAC',
                )
    tagElecAdv = (
                    "ENCUTGW","LASPH","LPEAD","LOPTICS","ENCUTGWSOFT"
                    'NKRED','NKREDX','NKREDY','NKREDZ',
                    'NGX','NGY','NGZ','NGXF','NGYF','NGZF',
                    'ODDONLY','EVENONLY','LASPH',
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
                    'LMIXTAU','NWRITE',
                    'AMIX','BMIX','TIME','SMASS',
                  )
    __tagNotImple= ()
    tagAll = tagElecBasic + tagElecAdv + tagElecXc + \
        tagIonBasic + tagIonSlab + tagPara + tagNotCateg
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

    def parse_tags(self, **keyval):
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

    def write(self, pathIncar="INCAR"):
        '''Write to INCAR file at 
        '''
        pass

    @classmethod
    def read_from_file(cls, pathIncar="INCAR"):
        '''Generate ``incar`` instance from INCAR file
        '''
        __kw = {}
        # TODO
        return cls(**__kw)

    @classmethod
    def minimal_electric(cls, pathIncar=None):
        pass

    @classmethod
    def minimal_ion_opt(cls, pathIncar=None):
        pass