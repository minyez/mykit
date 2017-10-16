#!/usr/bin/env python

from argparse import ArgumentParser
import subprocess as sp

description = '''
Generate structure files for parameter scanning and equation of state calculations
Use script optimize in wien2k. See details in wien2k userguide.

Currently support only fixed a/b/c ratio, i.e. varying volumes.
'''

parser = ArgumentParser(description=description)
parser.add_argument("-t", dest="stype", help="type of optimization",type=int, default=1)
parser.add_argument("-n", dest="nstruct", help="number of structures to calculate",type=int, default=11)
parser.add_argument("-s", dest="st", help="starting value. Integer",type=int, default=-10)
parser.add_argument("-e", dest="ed", help="ending value. Integer",type=int, default=10)
parser.add_argument("-D", dest="debug", help="debug mode",action='store_true')
opts = parser.parse_args()

ns = opts.nstruct
st = float(opts.st)
ed = float(opts.ed)
# in case of wrong input
if st > ed:
   st,ed = ed,st

dv = (ed-st)/float(ns-1)

ofile = 'opt.input'

with open(ofile,'w') as f:
    f.write(str(opts.stype)+'\n')
    f.write(str(ns)+'\n')
    for i in xrange(ns):
        f.write(str(round(st+float(i)*dv,1))+'\n')

sp.check_output('x optimize < %s' % ofile,stderr=sp.STDOUT,shell=True)
