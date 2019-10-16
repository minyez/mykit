#!/usr/bin/env python3
'''Generate structure files for parameter scanning and equation of state calculations

Use script optimize in wien2k. See details in wien2k userguide.
Currently only support 1D case, i.e. optimization type = 1,2,3,4 

    [1]  VARY VOLUME with CONSTANT RATIO A:B:C
    [2]  VARY C/A RATIO with CONSTANT VOLUME (tetr and hex lattices)
    [3]  VARY C/A RATIO with CONSTANT VOLUME and B/A (orthorh lattice)
    [4]  VARY B/A RATIO with CONSTANT VOLUME and C/A (orthorh lattice)
'''

import os
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import subprocess as sp
from shutil import copy2, rmtree
from mykit.wien2k.utils import get_casename
from mykit.core.utils import get_arith_prog


def pw_init_optimize_job():
    
    parser = ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-c", dest="casename", help="casename", default=None)
    parser.add_argument("-t", dest="stype", help="type of optimization.", type=int, default=1)
    parser.add_argument("-n", dest="nstruct", type=int, default=11, \
            help="number of structures to calculate")
    parser.add_argument("-s", dest="st", type=int, default=-10, \
            help="starting value. Integer")
    parser.add_argument("-e", dest="ed", type=int, default=10, \
            help="ending value. Integer")
    parser.add_argument("--ccmd", dest="f_calccmd", default=None, \
            help="file containing calculation command. Default: run_lapw -ec 0.0000001")
    parser.add_argument("--save", dest="savelapwdname", default=None, \
            help="Naming of the directory of savelapw command. Mustn't have space")
    parser.add_argument("--dir", dest="savedir", action="store_true", \
            help="also save each struct file to directoy V_x.xx/casename")
    # parser.add_argument("--run", dest="f_run", action='store_true', \
    #         help="flag for running the optimization")
    parser.add_argument("-D", dest="debug", action='store_true', \
            help="debug mode")
    
    args = parser.parse_args()

    if args.casename is None:
        cn = get_casename() 
    else:
        cn = args.casename
    
    ns = args.nstruct
    st = float(args.st)
    ed = float(args.ed)
    # in case of wrong input
    if st > ed:
        st, ed = ed, st
    
    if args.f_calccmd is None:
        calccmd = "run_lapw -ec 0.000001\n"
    else:
        with open(args.f_calccmd, 'r') as f_calccmd:
            calccmd = f_calccmd.readlines()

    vals = get_arith_prog(st, ed, n=ns)
    
    ofile = 'opt.input'
    
    with open(ofile, 'w') as f:
        f.write(str(args.stype)+'\n')
        f.write(str(ns)+'\n')
        for v in vals:
            f.write(str(round(v, 1))+'\n')
    
    sp.check_output('x optimize < %s' % ofile, stderr=sp.STDOUT, shell=True)
    
    # change MSR1 keyword to MSR1a in case.inm, in case of minimization of inner coordinates
    #try:
    #    sp.check_output("sed -i 's/MSR1 /MSR1a /g' %s.inm" % casename, \
    #       stderr=sp.STDOUT, shell=True)
    #except sp.CalledProcessError:
    #    print("WARNING: %s.inm is not found." % casename)
    
    # read the optimize.job
    with open('optimize.job', 'r') as h:
        lines_optjob = h.readlines()
    
    for i, line in enumerate(lines_optjob):
        words = line.split()
        len_line = len(words)
        if len_line == 0:
            continue
        elif words[0].startswith('#'):
            continue
        else:
            # write the calculation command from args.f_calccmd
            if words[0] == "run_lapw":
                lines_optjob[i] = ''.join(calccmd)
            # change the directory of savelapw
            elif words[0] == "save_lapw":
                if args.savelapwdname is not None:
                    lines_optjob[i] = "save_lapw -d ${i}_%s\n" % args.savelapwdname
                else:
                    lines_optjob[i] = "save_lapw -d ${i}_default\n"
    # back up
    copy2('optimize.job', 'optimize.job.old')
    
    with open('optimize.job', 'w') as h:
        for line in lines_optjob:
            h.write(line)
    
    # also save structures to directory
    if args.savedir:
        for i, v in enumerate(vals):
            dname = 'V_%4.2f' % (1+v*0.01)
            if os.path.isdir(dname):
                print("%s found. Remove old..." % dname)
                rmtree(dname)
            os.makedirs(os.path.join(dname, cn))
            fn = cn + '_' + {1: 'vol'}[args.stype] + f'{v:6.1f}'.replace(' ', '_') + '.struct'
            if os.path.isfile(fn):
                print("copying", fn, "...")
                copy2(fn, os.path.join(dname, cn, cn + '.struct'))


    # run the optimization
    # if args.f_run:
    #     fout = 'opt_type_%s.out' % args.stype
    #     ferr = 'opt_type_%s.error' % args.stype
    #     ofile = open(fout, 'w')
    #     efile = open(ferr, 'w')
    
    #     p = sp.Popen('./optimize.job', stdout=ofile, stderr=efile, shell=True)
    #     p.wait()
    
    #     ofile.close()
    #     efile.close()

if __name__ == "__main__":
    pw_init_optimize_job()
