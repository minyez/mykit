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

    __vaspTagval = {}

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

    def __init__(self, **incarargs):
        '''Initialize

        First, filter all incarargs to get all pwTags, xcTags, and remove their VASP equivalents,
        and then add those tags not implemented in the plane_wave tag mapping.
        '''
        kwargs = copy.deepcopy(incarargs)
        _vaspKs = []
        _vaspVs = []
        for _k in kwargs:
            if _k in self.__pwTagMaps:
                pass
            elif _k in self.tagAll:
                _v = kwargs.pop(_k)
        plane_wave.__init__(self, "vasp", **kwargs)

    def __parse_vasp_only_tags(self, **tags):
        pass

    def write(self, pathIncar="INCAR"):
        '''Write to INCAR file at 
        '''
        pass

    def add_tag(self, tagName, value):
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