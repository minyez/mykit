#!/usr/bin/env python
# coding=utf-8

# ====================================================
#
#     File Name : pw_classes.py
# Creation Date : 01-11-2017
# Last Modified : Wed 01 Nov 2017 03:11:39 PM CST
#    Created By : Min-Ye Zhang
#       Contact : stevezhang@pku.edu.cn
#       Purpose : classes for analysis of wien2k calculations
#
# ====================================================

import sys
from pw_anal_utils import Get_Casename

class w2k_struct():
    def __init__(self,casename=None):
        fn_struct = ''
        if casename is None:
            casename = Get_Casename()
        fn_struct = casename + '.struct'
        with open(fn_struct,'r') as f:
            self.strtlines = f.readlines()


# ==============================

