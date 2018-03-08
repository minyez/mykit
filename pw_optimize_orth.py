#!/usr/bin/env python
# coding=utf-8

# ====================================================
#
#     File Name : pw_optimize_orth.py
# Creation Date : 06-03-2018
#    Created By : Min-Ye Zhang
#       Contact : stevezhang@pku.edu.cn
#
# ====================================================

from __future__ import print_function
import sys, os, fnmatch
import subprocess as sp
from shutil import copy2
from argparse import ArgumentParser
from pw_anal_utils import Get_Casename
from pw_init_opt import pw_init_optimize_job
from pv_calc_utils import common_io_checkdir
from pv_anal_utils import vasp_anal_fit_EOS

def get_optimized_latt_const(casename, optdir, reference_struct):
    '''
    use equation of state to get the optimized boa/coa/vol
    casename: the casename
    optdir: boa, coa or vol
    reference_struct: the struct file to get the reference lattice constants
    '''
    # get the reference lattice constants
    h_ref = open(reference_struct,'r')
    for i, line in enumerate(h_ref):
        if i == 3:
            lattconst_ref_list = line.split()
        elif i >= 4:
            break
    h_ref.close()

    lattconst_ref_list = [float(x) for x in lattconst_ref_list]
    a_ref, b_ref, c_ref = lattconst_ref_list[0], lattconst_ref_list[1], \
                          lattconst_ref_list[2]
    boa_ref = b_ref/a_ref
    coa_ref = c_ref/a_ref
    vol_ref = a_ref * b_ref * c_ref

    # get the list of directories of saved calculations
    calcdir_list = []
    energy_list = []
    lattconst_list = []
    for idir in os.listdir('.'):
        if fnmatch.fnmatch(idir,casename+'_'+optdir+'_*'):
            if os.path.isdir(idir):
                calcdir_list.append(idir)

    item_list = []
    for x in calcdir_list:
        dir_suffix = x.split(casename+'_'+optdir)[1]
        for y in dir_suffix.split('_'):
            # append the first non-None element 
            if y != '':
                item_list.append(float(y))
                break
    
    #item_list = [float(x.split(casename+'_'+optdir)[1].split('_')[1]) for x in calcdir_list]

    for calcdir in calcdir_list:
        tote = 0
        os.chdir(calcdir)

        # get the total energy by finding :ENE
        with open(casename+'.scf','r') as h_scf:
            lines = h_scf.readlines()
        for line in lines:
            if line.startswith(":ENE"):
                tote = float(line.split()[-1])
        # save the last total energy
        energy_list.append(tote)

        # get the lattice constants
        h_struct = open(casename+'.struct','r')
        for i, line in enumerate(h_struct):
            if i == 3:
                lattconst_list.append([float(x) for x in line.split()])
            elif i >= 4:
                break
        h_struct.close()

        os.chdir('..')

    if optdir == 'boa':
        val_list = [x[1]/x[0] for x in lattconst_list]
    if optdir == 'coa':
        val_list = [x[2]/x[0] for x in lattconst_list]
    if optdir == 'vol':
        val_list = [x[0]*x[1]*x[2] for x in lattconst_list]

    # write the energy-val file.
    with open('../Ener_'+optdir,'w') as h_ener_val:
        h_ener_val.write("#ratio %s ene\n" % optdir)
        for i in xrange(len(val_list)):
            h_ener_val.write("%s %s %s\n" % (item_list[i],val_list[i],energy_list[i]))

    # get the optimized value of boa/coa/vol
    opte, optval, mod, modp = vasp_anal_fit_EOS('../Ener_'+optdir)
    if optdir == 'boa':
        a_opt = (vol_ref/coa_ref/optval)**(1./3.)
        b_opt = a_opt * optval
        c_opt = a_opt * coa_ref
    if optdir == 'coa':
        a_opt = (vol_ref/boa_ref/optval)**(1./3.)
        b_opt = a_opt * boa_ref
        c_opt = a_opt * optval
    if optdir == 'vol':
        optratio = (optval/vol_ref)**(1./3.)
        a_opt = a_ref * optratio
        b_opt = b_ref * optratio
        c_opt = c_ref * optratio

    print("    Optimized value: %10.6f" % optval)
    print("    Init. latt. const.: %10.6f%10.6f%10.6f" % (a_ref, b_ref, c_ref))
    print("    Optd. latt. const.: %10.6f%10.6f%10.6f" % (a_opt, b_opt, c_opt))
    return a_opt, b_opt, c_opt


def rewrite_latt_const(a, b, c, struct_old, struct_new):
    with open(struct_old, 'r') as h_s_old:
        lines = h_s_old.readlines()

    lines[3] = "%10.6f%10.6f%10.6f%10.6f%10.6f%10.6f\n" % (a,b,c,90,90,90)

    with open(struct_new, 'w') as h_s_new:
        for line in lines:
            h_s_new.write(line)


def perform_optimization(casename, optdir, init_lapw_cmd, options, starting_struct, optimized_struct):
    '''
    casename
    optdir: boa, coa or vol, corresponding to optimization of b/a, c/a and volume
    init_lapw_cmd: the batch command of init_lapw
    options: options for performing the function pw_init_optimize_job
    starting_struct: the path of the struct file to start
    optimized_struct: the path to save the optimized structure
    '''
    os.chdir(optdir)
    common_io_checkdir(casename)
    copy2('../'+starting_struct, casename+'/'+casename+'.struct')
    os.chdir(casename)

    # init_lapw
    sp.check_output(init_lapw_cmd, stderr=sp.STDOUT, shell=True)
    pw_init_optimize_job(options)

    # get the optimized lattice constants for the particular type of optimization
    a, b, c = get_optimized_latt_const(casename, optdir, '../../'+starting_struct)
    # rewrite the lattice constants in the struct file and save to a new one
    rewrite_latt_const(a, b, c, '../../'+starting_struct, '../../'+optimized_struct)

    os.chdir('../..')


def Main(ArgList):
    
    description='''
    Optimize an orthnormbic lattice with MSR1a method for minimization of inner coordinates.
    An optimization round includes: 1. optimize b/a with fixed volume and c/a;
    2. optimize c/a with fixed volume and optimized b/a;
    3. optimize volume with the optimized a:b:c ratio fixed.
    A couple of rounds of minimization should be executed to get a somewhat converged result.

    ATTENTION: a -8.0 eV is used as default.for core-valence splitting in init_lapw. Check if it is reasonable.
    Spin-polarized calculation is not supported yet.
    '''
    parser = ArgumentParser(description=description)
    parser.add_argument("-c", dest="casename", help="Case name", default=None)
    parser.add_argument("-n", dest="nrounds", help="The maximum number of optimization rounds",type=int, default=1)
    parser.add_argument("--vxc", dest="vxc", help="XC functional to use. Default: 13 (XC_PBE)", type=int, default=13)
    parser.add_argument("--nkp", dest="nkp", help="Total number of k-points in BZ", type=int, default=1000)
    parser.add_argument("--ecut", dest="ecut", help="Core-valence splitting in init_lapw. Default: -8.0 (Ry)", type=float, default=-8.0)
    parser.add_argument("--ccmd", dest="f_calccmd", help="file containing calculation command. Default: run_lapw -ec 0.000001", default=None)
    parser.add_argument("--save", dest="savelapwdname", help="Naming of the directory of savelapw command. Mustn't have space", default=None)
    parser.add_argument("-D", dest="debug",help="flag for debug mode.", action="store_true")

    opts = parser.parse_args()

    if opts.casename is None:
        casename = Get_Casename()
    else:
        casename = opts.casename

    case_struct            = casename + '.struct'
    case_round_init_struct = casename + '_init.struct'
    case_boa_optd_struct   = casename + '_boa_optd.struct'
    case_coa_optd_struct   = casename + '_coa_optd.struct'
    case_round_optd_struct = casename + '_optd.struct'

    if not os.path.exists(casename+'.struct'):
        print("ERROR: %s.struct is not found. Exit" % casename)
        sys.exit(1)

    # use a smaller number of structures (7) in b/a and c/a calculations
    # to avoid NN error (touching spheres)
    boa_options = ['-t','4','--run','-n','7','-s','-6','-e','6']
    coa_options = ['-t','3','--run','-n','7','-s','-6','-e','6']
    vol_options = ['-t','1','--run','-n','9','-s','-8','-e','8']

    refined_ArgList = ArgList
    # remove the options that pw_init_optimize_job does not have
    # explicitly, --vxc, --nkp, --ecut
    for opt in ['--vxc','--nkp','--ecut']:
        if opt in refined_ArgList:
            opt_index = refined_ArgList.index(opt)
            # delete the tag
            del refined_ArgList[opt_index]
            # delete the value
            del refined_ArgList[opt_index]

    # change the cal_ccmd file to the absolute path
    if '--ccmd' in refined_ArgList:
        ccmd_index = refined_ArgList.index('--ccmd')
        refined_ArgList[ccmd_index+1] = os.path.abspath(opts.f_calccmd)

    # set the batch init_lapw command
    init_lapw_cmd = 'init_lapw -b' + ' -ecut %s'%opts.ecut + ' -vxc %s'%opts.vxc + ' -numk %s'%opts.nkp

    try:
        assert opts.nrounds > 0
    except:
        raise ValueError("nrounds should be positive.")

    #return

    for i in xrange(opts.nrounds):
        round_dir = 'round-%s'%(i+1)
        round_dir_prev = 'round-%s'%i
        common_io_checkdir(round_dir)

        if i == 0:
            copy2(case_struct,round_dir+'/'+case_round_init_struct)
        else:
            copy2(round_dir_prev+'/'+case_round_optd_struct, \
                  round_dir+'/'+case_round_init_struct)

        os.chdir(round_dir)
        common_io_checkdir('boa')
        common_io_checkdir('coa')
        common_io_checkdir('vol')
        
        # perform b/a optimization
        print("Start optimization Round %d" % (i+1))
        print("  Type: boa")
        perform_optimization(casename, 'boa', init_lapw_cmd, refined_ArgList.extend(boa_options), \
                             case_round_init_struct, case_boa_optd_struct)

        # perform c/a optimization
        print("  Type: coa")
        perform_optimization(casename, 'coa', init_lapw_cmd, refined_ArgList.extend(coa_options), \
                             case_boa_optd_struct, case_coa_optd_struct)

        # perform EOS calculation with fixed abc ratio
        print("  Type: vol")
        perform_optimization(casename, 'vol', init_lapw_cmd, refined_ArgList.extend(vol_options), \
                             case_coa_optd_struct, case_round_optd_struct)

        # finish of one optimization round
        print("Finish optimization Round %d" % (i+1))
        os.chdir('..')


# ==============================

if __name__ == "__main__":
    Main(sys.argv)
