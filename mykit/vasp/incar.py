# coding = utf-8
'''this module implements the ``incar`` class that manages the VASP input file INCAR
'''
import copy
from mykit.core.planewave import plane_wave, planewaveError
from mykit.core.xc import xc_control, xcError

class incarError(Exception):
    pass


class incar(plane_wave, xc_control):
    '''manage tags and IO of VASP input file INCAR
    '''
    # For VASP tags that is not easily to find analogy in other programs
    vaspOnlyTags = ()
    tagElecBasic = (
                    "ISTART","LWAVE","IALGO","EDIFF","ENCUT","NBANDS",
                    "ICHARG","LCHARG","PREC","ALGO","IALGO","ISPIN","NELM",
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
    tagAll = tagElecBasic + tagElecAdv + tagElecXc + \
        tagIonBasic + tagIonSlab + tagPara + tagNotCateg
    tagMap2Base = {}
    __availMap2Pw = plane_wave.map2pwtags(*tagAll, progFrom="vasp")
    # TODO xc map
    __availMap2Xc = (None,)*len(tagAll)
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
        # _vaspKs = []
        # _vaspVs = []
        # for _k in kwargs:
        #     if _k in self.__pwTagMaps:
        #         pass
        #     elif _k in self.tagAll:
        #         _v = kwargs.pop(_k)
        super(incar, self).__init__(**incarArgs)
        self.__parse_vasptags(**incarArgs)

    def parse_tags(self, **keyval):
        super(incar,self).parse_tags(**keyval)
        self.__parse_vasptags(**keyval)

    def __parse_vasptags(self, **tags):
        __keyFiltered = []
        if len(tags) != 0:
            __keyFiltered = self.__filter_tags_incar_not_pw_xc(*tags.keys())
        for _k in __keyFiltered:
            self.__vaspTags.update({_k:tags[_k]})

    @classmethod
    def __filter_tags_incar_not_pw_xc(cls, *tags):
        '''Get tags that are implemented in incar but not implemented as common plane_wave or xc_control tags

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
        # Filter xcTags and pwTags
        for _k, _v in __tagsDict.items():
            if _k in cls.tagMap2Base.values():
                pass
            elif _k in cls.tagAll:
                _filter1.append(_k)
        if len(_filter1) == 0:
            return _filter1
        # Filter tags which can be mapped to xcTags and pwTags
        # __map2pw = plane_wave.map2pwtags(*_filter1, progFrom="vasp")
        # __map2xc = [None,] * len(_filter1)
        # # __map2xc = xc_control.map2xctags(*_filter1, progFrom="vasp")
        _filter2 = []
        for _t in _filter1:
            if _t not in cls.tagMap2Base:
                _filter2.append(_t)
        return _filter2

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