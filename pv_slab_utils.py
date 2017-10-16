#!/usr/bin/env python
#
# ver.0.0.1 by M.-Y. Zhang
# ====================================================
# ====================================================

import sys, os, shutil, subprocess, copy
import commands, string
from argparse import ArgumentParser

# ====================================================

def vasp_get_NPAR(np):

    npar = 1
    nsqrt_max = int(sqrt(np)) + 1
    nsqrt = nsqrt_max
    while nsqrt > 0:
        if nsqrt % 2 == 0 and np % nsqrt == 0:
            npar = nsqrt
            break
        nsqrt = nsqrt - 1
    return npar
