#!/usr/bin/env python

from argparse import ArgumentParser
from pw_utils import w2k_get_casename
import subprocess as sp
import os,sys

description =  '''
Read the observable-volume data from an optimize.job execution.
SCF data are supposed to save in a directory named as claimed in the optimize.job script, i.e. save_lapw.
'''

obs_list = ['ENE','GAP','VOL','FER',\
            'LAT1','LAT2','LAT3']
loc      = [    9,    6,    7,   10,\
                5 ,    6 ,    7 ]

parser = ArgumentParser(description=description)
parser.add_argument("-o", dest="obs", help="the observable to monitor. Available: %s" % ','.join(obs_list),default="ENE")
parser.add_argument("-D", dest="debug", help="debug mode",action='store_true')
parser.add_argument("-p", dest="save", help="plot data (fitting BM-EOS)",action='store_true')
opts = parser.parse_args()
obs = opts.obs

if not obs in obs_list:
    print "Unsupported observable monitor. Exit"
    sys.exit(1)

ifile = 'optimize.job'
ofile = 'Ener_Vol'
struct_flag = False
obj_struct = [ ]
casename = w2k_get_casename()
if opts.debug: print casename


with open(ifile,'r') as f:
    lines = f.readlines()
    for line in lines:
        words = line.split()
        if len(words) == 0:
            continue
        if words[0] == ')':
            struct_flag = False
        if struct_flag:
            obj_struct.append(words[0])
        if words[0] == 'foreach':
            struct_flag = True
        if words[0] == 'save_lapw':
            dirname = words[words.index('-d')+1]


for x in xrange(len(obj_struct)):
    obj_struct[x] = obj_struct[x] + dirname[4:]

if opts.debug: print obj_struct

vol_data = []
obs_data = []

for struct in obj_struct:
    try:
        os.chdir(struct)
    except:
        continue
    obs1 = sp.check_output("awk '/:%s/ {print $%i}' %s.scf | tail -1" % (obs,loc[obs_list.index(obs)],casename),shell=True).strip()
    vol1 = sp.check_output("awk '/:%s/ {print $%i}' %s.scf | tail -1" % ('VOL',loc[obs_list.index('VOL')],casename),shell=True).strip()
    vol_data.append(vol1)
    obs_data.append(obs1)
    if opts.debug: print vol1,obs1
    os.chdir('..')


with open(ofile,'w') as fout:
    fout.write("# Vol %s\n" %obs)
    for i in xrange(len(vol_data)):
        fout.write("%3i   %s   %s\n" % (i+1,vol_data[i],obs_data[i]))
