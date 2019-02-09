#!/usr/bin/env python
# generate a number of structures with different volumes for plotting energy-volume curve

from argparse import ArgumentParser
import numpy as np
import os
import errno

def create_path(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

parser = ArgumentParser()
parser.add_argument("-i", help="input struct file, default POSCAR")
parser.add_argument("-s", help="start volume ratio, default 0.92", type=float)
parser.add_argument("-e", help="end volume ratio, default 1.08", type=float)
parser.add_argument("-d", help="interval of volume ratio, default 0.02", type=float)
args = parser.parse_args()

if args.i:
    struct_file = args.i
else:
    struct_file = "POSCAR"

if args.s:
    rvs = args.s
else:
    rvs = 0.92

if args.e:
    rve = args.e
else:
    rve = 1.08

if args.d:
    interval = args.d
else:
    interval = 0.02

with open(struct_file,'r') as f:
    lines = f.readlines()

scale_old = float(lines[1].split()[0])
grid = round(( rve - rvs ) / interval) + 1
struct_range_np = np.linspace(rvs,rve,grid)
struct_range_str = [ ]
for x in struct_range_np:
     struct_range_str.append("%3.2f" % x)
#print rvs,rve,interval,(rve-rvs)/interval,round((rve-rvs)/interval),grid
#print struct_range_np
#print struct_range_str

# create directories with POSCAR for each volume
for x in struct_range_str:
    create_path("./V_%s" % x)
    ra = np.power(float(x),1/3.0)
    lines[1]= "%s\n" % str(ra * scale_old)
    os.chdir("./V_%s" % x)
    with open("POSCAR",'w') as poscar:
        for line in lines:
            poscar.write(line)
    os.chdir("../")

# write bash code to run the scanning
#with open("run-Ener-Vol.sh",'w') as f:
#    f.write("#!/usr/bin/env bash")
#    f.write("# vasp=")
#    f.write("for i in %s" % " ".join())
#    f.write("do")
#    f.write()
#    f.write("done")


