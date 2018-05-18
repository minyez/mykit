#!/usr/bin/env python

from __future__ import print_function
from argparse import ArgumentParser
import subprocess as sp
import sys
from shutil import copy2
from pw_anal_utils import w2k_get_casename


def pw_init_optimize_job(ArgList):
    description = '''
    Generate structure files for parameter scanning and equation of state calculations
    Use script optimize in wien2k. See details in wien2k userguide.
    Currently only support 1D case, i.e. optimization type = 1,2,3,4 
        [1]  VARY VOLUME with CONSTANT RATIO A:B:C
        [2]  VARY C/A RATIO with CONSTANT VOLUME (tetr and hex lattices)
        [3]  VARY C/A RATIO with CONSTANT VOLUME and B/A (orthorh lattice)
        [4]  VARY B/A RATIO with CONSTANT VOLUME and C/A (orthorh lattice)
    '''
    
    parser = ArgumentParser(description=description)
    parser.add_argument("-c", dest="casename", help="casename", default=None)
    parser.add_argument("-t", dest="stype", help="type of optimization.",type=int, default=1)
    parser.add_argument("-n", dest="nstruct", help="number of structures to calculate",type=int, default=11)
    parser.add_argument("-s", dest="st", help="starting value. Integer",type=int, default=-10)
    parser.add_argument("-e", dest="ed", help="ending value. Integer",type=int, default=10)
    parser.add_argument("--ccmd", dest="f_calccmd", help="file containing calculation command. Default: run_lapw -ec 0.000001", default=None)
    parser.add_argument("--save", dest="savelapwdname", help="Naming of the directory of savelapw command. Mustn't have space", default=None)
    parser.add_argument("--run", dest="f_run", help="flag for running the optimization",action='store_true')
    parser.add_argument("-D", dest="debug", help="debug mode",action='store_true')
    
    opts = parser.parse_args(ArgList[1:])
   
    if opts.debug:
        print(ArgList)

    if opts.casename == None:
        casename = w2k_get_casename() 
    else:
        casename = opts.casename
    
    ns = opts.nstruct
    st = float(opts.st)
    ed = float(opts.ed)
    # in case of wrong input
    if st > ed:
       st,ed = ed,st
    
    if opts.f_calccmd == None:
        calccmd = "run_lapw -ec 0.000001\n"
    else:
        with open(opts.f_calccmd,'r') as f_calccmd:
            calccmd = f_calccmd.readlines()
    if ns != 1:
        dv = (ed-st)/float(ns-1)
    else:
        dv = 0.0
    
    ofile = 'opt.input'
    
    with open(ofile,'w') as f:
        f.write(str(opts.stype)+'\n')
        f.write(str(ns)+'\n')
        for i in xrange(ns):
            f.write(str(round(st+float(i)*dv,1))+'\n')
    
    sp.check_output('x optimize < %s' % ofile,stderr=sp.STDOUT,shell=True)
    
    # change MSR1 keyword to MSR1a in case.inm, in case of minimization of inner coordinates
    #try:
    #    sp.check_output("sed -i 's/MSR1 /MSR1a /g' %s.inm" % casename, stderr=sp.STDOUT, shell=True)
    #except sp.CalledProcessError:
    #    print("WARNING: %s.inm is not found." % casename)
    
    # read the optimize.job
    with open('optimize.job','r') as h_optimize:
        lines_optjob = h_optimize.readlines()
    
    for i in xrange(len(lines_optjob)):
        line = lines_optjob[i]
        words = line.split()
        len_line = len(words)
        if len_line == 0:
            continue
        elif words[0].startswith('#'):
            continue
        else:
            # write the calculation command from opts.f_calccmd
            if words[0] == "run_lapw":
                lines_optjob[i] = ''.join(calccmd)
            # change the directory of savelapw
            elif words[0] == "save_lapw":
                if opts.savelapwdname is not None:
                    lines_optjob[i] = "save_lapw -d ${i}_%s\n" % opts.savelapwdname
                else:
                    lines_optjob[i] = "save_lapw -d ${i}_default\n"
    
    copy2('optimize.job','optimize.job.old')
    
    with open('optimize.job','w') as h_optimize_new:
        for line in lines_optjob:
            h_optimize_new.write(line)

    # run the optimization
    if opts.f_run:
        fout = 'opt_type_%s.out' % opts.stype
        ferr = 'opt_type_%s.error' % opts.stype
        ofile = open(fout,'w')
        efile = open(ferr,'w')
    
        p=sp.Popen('./optimize.job',stdout=ofile,stderr=efile,shell=True)
        p.wait()
    
        ofile.close()
        efile.close()


if __name__ == "__main__":
    pw_init_optimize_job(sys.argv)
