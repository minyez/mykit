#!/usr/bin/env python
# coding=utf-8

# ====================================================
#
#     File Name : pw_classes.py
# Creation Date : 01-11-2017
# Last Modified : Tue 07 Nov 2017 07:25:37 PM CST
#    Created By : Min-Ye Zhang
#       Contact : stevezhang@pku.edu.cn
#       Purpose : classes for analysis of wien2k calculations
#
# ====================================================

import sys
from pw_anal_utils import w2k_get_casename
#from w2k_utils import w2k_get

class w2k_struct():
    def __init__(self,casename=None):
        fn_struct = ''
        if casename is None:
            casename = w2k_get_casename()
        self.casename = casename
        fn_struct = casename + '.struct'
        with open(fn_struct,'r') as f:
            self.strtlines = f.readlines()


# ==============================

class w2k_energy_band():
    def __init__(self,ifile=None,nat=None):
        if ifile is None:
            self.ifile = w2k_get_casename()+'.energy-band'
        else:
            self.ifile = ifile
        if nat is None:
            pass
            #self.nat=w2k_get(casename,'nat')
        else:
            self.nat = nat
