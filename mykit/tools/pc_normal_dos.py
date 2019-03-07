#!/usr/bin/env python
# coding=utf-8

# ====================================================
#
#     File Name : pc_normal_dos.py
# Creation Date : 07-11-2017
# Last Modified : Wed 08 Nov 2017 08:25:39 PM CST
#    Created By : Min-Ye Zhang
#       Contact : stevezhang@pku.edu.cn
#       Purpose : Common utility to normalize with respect to the highest peak,
#                 also possible to shift the dos.
#
# ====================================================

from argparse import ArgumentParser
import os,sys

def Main(ArgList):
# =================== Parser ==========================
    norm_suffix = '_norm'

    description = '''
    Common utility to normalize with respect to the highest peak and shift the dos.
    '''
    parser = ArgumentParser(description=description)
    parser.add_argument('ifiles',nargs='+',help="data files")
    parser.add_argument("-s",dest='shifts',nargs='+',help="shift of data file",type=float)
    parser.add_argument("-r",dest='range',nargs=2,help="range of shift",type=float,default=[-1000,1000])
    parser.add_argument("-c",dest='calibrate',nargs=2,help="two points to linearly calibrate the baseline",type=float)
    parser.add_argument("-D",dest='debug',help="flag for debug mode",action='store_true')

    opts = parser.parse_args()
    ifiles = opts.ifiles
    shifts = opts.shifts
    peak_locs = [0.0]*len(ifiles)
    baseline = [0.0]*len(ifiles)
    peakval = [0.0]*len(ifiles)
    peaks = []
    if opts.calibrate is not None:
        cbp1 = min(opts.calibrate)
        cbp2 = max(opts.calibrate)

    if shifts is None:
        shifts = [0.0]*len(ifiles)
    elif len(ifiles) > len(shifts):
        shifts.extend([0.0]*(len(ifiles)-len(shifts)))
    elif len(ifiles) < len(shifts):
        shifts = shifts[0:len(ifiles)]
#    print ifiles
#    print shifts

    for i in xrange(len(ifiles)):
        peak_temp = []
        ifilename = ifiles[i]
        ener = []
        dos  = []
        maxval = -1.0
        minval = 10000
        loc_norm_peak = 0.0
        with open(ifilename,'r') as if_dos:
            dos_lines = if_dos.readlines()
        for line in dos_lines:
            line_str = line.strip().split()
            if len(line_str) == 0:
                continue
            if line_str[0].startswith('#'):
                continue
            ener1 = float(line_str[0])
            dos1 = float(line_str[1])
            if line_str[1] == 'NaN':
                continue
            if shifts[i] != 0.0 and ((ener1-opts.range[0])*(ener1-opts.range[1])<=0.0):
                ener.append(ener1+shifts[i])
            else:
                ener.append(ener1)
            dos.append(dos1)
            if dos1 > maxval:
                maxval = dos1
                loc_norm_peak = ener1
            if dos1 < minval:
                minval = dos1
# calibrate the base line and normalize the data
        dos = [(x-minval)/(maxval-minval)*100.0 for x in dos]
        baseline[i] = minval
        peakval[i] = maxval
        peak_locs[i] = loc_norm_peak

    for i in xrange(len(ifiles)):
        ifilename = ifiles[i]
        if not ifilename.endswith(norm_suffix):
            ofilename = ifilename + norm_suffix
        else:
            ofilename = ifilename
        with open(ofilename,'w') as of_dos:
            of_dos.write("# original peak at %10.5f\n" % peak_locs[i])
            of_dos.write("# original peak and baseline: %10.5f,%10.5f\n" % (peakval[i],baseline[i]))
#            print shifts[i]
            if shifts[i] != 0.0:
                of_dos.write("# shifted with %10.5f (range %6.2f to %6.2f)\n" % (shifts[i],min(opts.range),max(opts.range)))
            else:
                of_dos.write("# unshifted\n")
            for i in xrange(len(dos)):
                of_dos.write("%10.5f %15.7f\n" % (ener[i],dos[i]))



# ==============================

if __name__ == "__main__":
    Main(sys.argv)

