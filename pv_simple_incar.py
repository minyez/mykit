#!/usr/bin/env python

from pv_calc_utils import *

with open('INCAR','w') as incar:
    vasp_write_incar_minimal_elec(incar,tag_xc="PBE")
