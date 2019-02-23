# coding = utf-8
'''this module implements the ``incar`` class that manages the VASP input file INCAR
'''

from mykit.core.planewave import plane_wave


class incar(plane_wave):

    __tagElecBasic = {
        "ISTART" : ((int), 0),
        "LWAVE" : ((bool), True)
    }
    __tagElecAdv = {}
    __tagIonBasic = {}
    __tagIonSlab = {}
    # ['LCHARG','LWAVE','LDIPOL','LVHAR','LVTOT','LHFCALC','LASPH','LMIXTAU']
    # ['NGX','NGY','NGZ','NGXF','NGYF','NGZF','NBANDS','ISTART',\
    #   'ICHARG','ISPIN','NWRITE','NELM','NSW','IBRION','ISIF','ISYM',\
    #   'IALGO','ISMEAR','NKRED','NKREDX','NKREDY','NKREDZ','ODDONLY',\
    #   'EVENONLY','LMAXTAU','IDIPOL','NELMIN','NPAR']
    # ['AEXX','AGGAX','AGGAC','ALDAC','ENCUT','ENCUTGW','EDIFF',\
    #   'EDIFFG','POTIM','NELECT','SMASS','AMIX','BMIX','SIGMA',\
    #   'HFSCREEN','TIME']
    # ['PREC','PRECFOCK','ALGO','LREAL','GGA','METAGGA']
    # ['MAGMOM','DIPOL']

    xc_dict = {
                # 'LDA'    : {'GGA'     : 'CA'} ,
                # 'PBE'    : {'GGA'     : 'PE'} ,
                # 'PBESOL' : {'GGA'     : 'PS'} ,
                # 'RPBE'   : {'GGA'     : 'RP'} ,
                # 'SCAN'   : {'METAGGA' : 'SCAN', 'LASPH' : '.TRUE.'} ,
                # 'PBE0'   : self.tag_PBE0  ,
                # 'HSE06'  : self.tag_HSE06 ,
                # 'HF'     : self.tag_HF
                }

    def __init__(self, **kwargs):

        super(incar, self).__init__("vasp", **kwargs)

    @classmethod
    def read_from_file(cls, pathIncar):
        '''Generate ``incar`` instance from
        '''
        __kw = {}
        return cls(**__kw)