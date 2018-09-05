#!/usr/bin/env python

from __future__ import print_function
import sys, os
from shutil import copy2
from argparse import ArgumentParser
from pv_calc_utils import vasp_write_incar_minimal_elec, vasp_io_change_tag, \
                          vasp_write_kpoints_basic, vasp_write_incar_exact, \
                          vasp_vasprun_zmy, vasp_vaspcmd_zmy, vasp_io_get_NPAR
from pv_anal_utils import vasp_anal_get_outcar
from pv_classes import vasp_read_poscar
from pc_utils import common_io_checkdir, common_io_cleandir

# =====================================================

# Step 1: self-consistent calculation of charge density
def step_1_scf(chg_dir, tag_xc, encut, nks,\
               npar, lmm, vasp_cmd, f_metal=False, smear=[0, 0.05]):
    if not os.path.exists(chg_dir):
        common_io_checkdir(chg_dir)
    else:
        print("Charge calculation done before. Restart")
        common_io_cleandir(chg_dir)
    copy2('POSCAR', chg_dir)
    copy2('POTCAR', chg_dir)
    os.chdir(chg_dir)

    with open('INCAR', 'w') as incar:
        vasp_write_incar_minimal_elec(incar, tag_xc, encut, npar=npar, mode_smear=smear)
    if f_metal:
        vasp_io_change_tag('INCAR', 'SIGMA', new_val=0.2)
    vasp_io_change_tag('INCAR', 'LMAXMIX', new_val=lmm)

    vasp_write_kpoints_basic(nks)
    vasp_vasprun_zmy(vasp_cmd,'out','error')
    # get the maximum number of planewaves
    mnpw = vasp_anal_get_outcar('mnpw')
    os.chdir('..')
    return mnpw

# ==================================================

# Step 2: non-self-consistent calculation of unoccupied orbitals
def step_2_exact(chg_dir, exact_dir, tag_xc, encut, nks, nks_gw, nbands, nedos, lwannier,\
                 npar, lmm, vasp_cmd, f_metal=False, smear=[0, 0.05]):
    common_io_cleandir(exact_dir)
    copy2('POSCAR', exact_dir)
    copy2('POTCAR', exact_dir)
    os.chdir(exact_dir)
    vasp_write_kpoints_basic(nks_gw)
    with open('INCAR', 'w') as incar:
        vasp_write_incar_exact(incar, tag_xc, encut, nb=nbands, npar=npar, mode_smear=smear, lwannier=lwannier)
    vasp_io_change_tag('INCAR', 'NEDOS', new_val=nedos)
    if f_metal:
        vasp_io_change_tag('INCAR', 'SIGMA', new_val=0.2)

# if kpoints for the first and second calculation is the same, then it is not necessary to fix charge density, and the wavefunction in Step 1 can be used
    if not nks == nks_gw:
        copy2('../'+chg_dir+'/CHGCAR', '.')
        vasp_io_change_tag('INCAR', 'ISTART')
        vasp_io_change_tag('INCAR', 'ICHARG', new_val=11)
        vasp_io_change_tag('INCAR', 'LMAXMIX', new_val=lmm)
    else:
        copy2('../'+chg_dir+'/WAVECAR','.')
        vasp_io_change_tag('INCAR','ICHARG')
        vasp_io_change_tag('INCAR','ISTART',new_val=1)

    vasp_vasprun_zmy(vasp_cmd,'out','error')
    os.chdir('..')

# ==================================================

# Step 3: GW calculation
def step_3_gw(exact_dir,gw_dir,tag_xc,encut, encutgw, nbands,nedos,nomega,lwannier,\
              vasp_cmd,gw_mode="G0W0",f_metal=False,smear=[0,0.05]):
# need to add more rules
    common_io_cleandir(gw_dir)
    copy2('POSCAR',gw_dir)
    copy2('POTCAR',gw_dir)
    copy2(exact_dir+'/KPOINTS',gw_dir)
    copy2(exact_dir+'/WAVECAR',gw_dir)
    if not f_metal:
        copy2(exact_dir+'/WAVEDER',gw_dir)
    os.chdir(gw_dir)

    with open('INCAR','w') as incar:
        vasp_write_incar_exact(incar,tag_xc,encut,nb=nbands,npar=1,mode_smear=smear,lwannier=lwannier)
        incar.write(' NEDOS = %s\n' % nedos )
        incar.write(' NOMEGA = %s\n' % nomega )
    if encutgw != 0:
        vasp_io_change_tag('INCAR','ENCUTGW',new_val=encutgw)
    if gw_mode == "G0W0":
        vasp_io_change_tag('INCAR','ALGO',new_val='GW0')
    if f_metal:
        vasp_io_change_tag('INCAR','SIGMA',new_val=0.2)

    vasp_vasprun_zmy(vasp_cmd,'out','error')
    os.chdir('..')

# =====================================================

def Main(ArgList):

    description = '''
    Prepare INCAR files and run the GW calculation
    Need to first prepare POSCAR and POTCAR
    '''
    testmode_avail = ['ENCUTGW','NBANDS']

    chg_dir = '1_charge'
    exact_dir = '2_exact'
    gw_dir = '3_gw'
    conv_egw_dir = 'CONV_encutgw'
    conv_nbs_dir = 'CONV_nbands'

# =================== Parser ==========================

    parser = ArgumentParser(description=description)
    parser.add_argument("-e",dest='encut',type=float,default=0,help="Planewave cutoff. 0 will set largest ENMAX")
    parser.add_argument("--egw",dest='encutgw',type=float,default=0,help="dielectric cutoff for GW calculation, i.e. ENCUTGW")
    parser.add_argument("-n",dest='nproc',type=int,default=1,help="Number of processors ")
    parser.add_argument("-b",dest='nbands',default=None,help="Number of bands")
    parser.add_argument("-x",dest='tag_xc',default=None,help="type of XC functional for input orbitals, \
                        None for LEXCH in POTCAR")
    parser.add_argument("-k",dest='nk',default=[6,6,6],help="Number of kpoints for SCF, format like [kx,ky,kz]")
    parser.add_argument("--kg",dest='nk_gw',default=[6,6,6],help="Number of kpoints for GW, format like [kgx,kgy,kgz]")
    parser.add_argument("-m",dest='f_metal',action = 'store_true',help="flag for metal system")
    parser.add_argument("-w",dest='nomega',default=50,help="Number of frequency points. Default 50")
    parser.add_argument("-D",dest='debug',action='store_true',help="flag for debug mode")
    parser.add_argument("-v",dest='vasp_path',default="vasp",help="Path of vasp executive")
    parser.add_argument("--gw",dest='gw_mode',default="G0W0",help="Self-consistent level of GW. Default G0W0")
    parser.add_argument("--restart",dest='chg_start',action='store_true',help="flag for restarting from SCF")
    parser.add_argument("--wannier",dest='lwannier',action='store_true',help="flag for wannier90 to plot band structure")
    parser.add_argument("--test",dest='testmode',default=None,help="flag for convergence test with respect to testmode,\
                        e.g. ENCUTGW, NBANDS")

    opts = parser.parse_args()
    npar = vasp_io_get_NPAR(opts.nproc)
    nks = opts.nk
    nks_gw = opts.nk_gw
    if opts.debug:
        print(nks)
        print(nks_gw)
        sys.exit(0)
    if opts.nbands is not None:
        nbands = opts.nbands
    else:
# setting default, i.e. 8 times the number of atoms times nproc
# DEPRECATED
#        poscar = vasp_read_poscar()
#        nbands = 8*poscar.natoms*opts.nproc
        nbands = 0

    testmode = None
    if opts.testmode is not None:
        if opts.testmode.upper() in testmode_avail:
            testmode = opts.testmode.upper()
        else:
            print("%s is not supported now. Test mode off." % opts.testmode)

# grid of density of states
    nedos = 3000
# lmaxmix, change to larger value if f-element present
    lmm = 4
# set vasp variable
    vasp_path, vasp_cmd = vasp_vaspcmd_zmy(opts.nproc, 'mpirun', opts.vasp_path)

# ================== End Parser =========================

# 1) run an SCF-calculation for the charge density
    print("== Step 1: SCF for charge density ==")
#     if CHGCAR does not exist but want to start from CHGCAR, report error
    if not os.path.exists(chg_dir) and opts.chg_start:
        print("Charge has not yet been calculated. Exit.")
        sys.exit(1)
#     if CHGCAR calculated and you do not want to start from SCF
    elif os.path.exists(chg_dir+'/CHGCAR') and not opts.chg_start:
        print("Charge calculation done. Move on")
        print("(Hint: Better check past setting of ENCUT and k-points)")
#     CHGCAR do not exist, or want to redo CHGCAR calculation
    else:
        mnpw = step_1_scf(chg_dir, tag_xc=opts.tag_xc, encut=opts.encut, nks=nks, npar=npar,\
                  vasp_cmd=vasp_cmd, f_metal=opts.f_metal, lmm=lmm, smear=[0, 0.05])
        if nbands == 0:
            nbands = mnpw

# IF testmode is not set
    if not testmode:
# 2) run a non-SCF calculation to get desired number of bands
# 3) run a GW calculation with WAVECAR and WAVEDER from Step 2
        print("== Step 2: Non-SCF for unc. bands ==")
        step_2_exact(chg_dir,exact_dir,opts.tag_xc,opts.encut,nks,nks_gw,nbands,nedos,lwannier=opts.lwannier,\
                 npar=npar,lmm=lmm,vasp_cmd=vasp_cmd,f_metal=opts.f_metal,smear=[0,0.05])

        print("== Step 3: GW                     ==")
        step_3_gw(exact_dir,gw_dir,opts.tag_xc,opts.encut, opts.encutgw, nbands=nbands,nedos=nedos,nomega=opts.nomega,\
                  lwannier=opts.lwannier,vasp_cmd=vasp_cmd,gw_mode=opts.gw_mode,f_metal=opts.f_metal,smear=[0,0.05])

# IF testmode is set to ENCUTGW
    elif testmode == 'NBANDS':
        common_io_cleandir(conv_nbs_dir)
        mnpw = int(vasp_anal_get_outcar('mnpw',outcar=chg_dir+'/OUTCAR'))
        nbands_scf = int(vasp_anal_get_outcar('nb',outcar=chg_dir+'/OUTCAR'))
        copy2('POSCAR',conv_nbs_dir)
        copy2('POTCAR',conv_nbs_dir)
        os.chdir(conv_nbs_dir)
        nbands_list = []
        for x in xrange(2, 8):
            nbands_ori = nbands_scf + x*int((mnpw-nbands_scf)/10.0)
            nbands_ori = int(nbands_ori/opts.nproc+1)*opts.nproc
            if nbands_ori not in nbands_list and (nbands_ori<mnpw):
                nbands_list.append(nbands_ori)

        for nb in nbands_list:
            print("== Step 2: Non-SCF for unc. bands ==")
            step_2_exact('../'+chg_dir,exact_dir,opts.tag_xc,opts.encut,nks,nks_gw,nbands=nb,nedos=nedos,\
                    lwannier=opts.lwannier, npar=npar,lmm=lmm,vasp_cmd=vasp_cmd,f_metal=opts.f_metal,smear=[0,0.05])
            print("== Step 3:    GW   BANDS:%5s    ==" % nb)
            step_3_gw(exact_dir,gw_dir+'_nb_'+str(nb),opts.tag_xc,opts.encut,nbands=nb,nedos=nedos,nomega=opts.nomega,\
                  lwannier=opts.lwannier,vasp_cmd=vasp_cmd,gw_mode=opts.gw_mode,f_metal=opts.f_metal,smear=[0,0.05])

        gap_list = []
        for nb in nbands_list:
            os.chdir(gw_dir+'_nb_'+str(int(nb)))
            gap_list.append(vasp_anal_get_outcar('gap'))
            os.chdir('..')
        os.chdir('..')

        #write to output file
        if os.path.exists('Gap_NBANDS.dat'):
            copy2('Gap_NBANDS.dat','Gap_NBANDS.dat_old')
        with open('Gap_NBANDS.dat', 'w') as f:
            f.write("#MNPW=%6s\n" % mnpw)
            f.write("#NBANDS gap\n")
            for i in xrange(len(nbands_list)):
                f.write("%5.1f %12.6f\n" % (nbands_list[i], gap_list[i]))

# IF testmode is set to ENCUTGW
    elif testmode == 'ENCUTGW':
        common_io_cleandir(conv_egw_dir)
        mnpw = int(vasp_anal_get_outcar('mnpw', outcar=chg_dir+'/OUTCAR'))
        nbands_scf = int(vasp_anal_get_outcar('nb', outcar=chg_dir+'/OUTCAR'))
        encut = int(vasp_anal_get_outcar('encut', outcar=chg_dir+'/OUTCAR'))

        if opts.nbands is not None:
            nbands = opts.nbands
        else:
            nbands = int(mnpw/2/opts.nproc+1)*opts.nproc

        copy2('POSCAR',conv_egw_dir)
        copy2('POTCAR',conv_egw_dir)
        os.chdir(conv_egw_dir)
        print("== Step 2: Non-SCF for unc. bands ==")
        step_2_exact('../'+chg_dir,exact_dir,opts.tag_xc,opts.encut,nks,nks_gw,nbands=nbands,nedos=nedos,lwannier=opts.lwannier,\
                 npar=npar,lmm=lmm,vasp_cmd=vasp_cmd,f_metal=opts.f_metal,smear=[0,0.05])
        egw_list = [encut/10.0*x for x in xrange(2,7)] # up to half of the ENCUT
        for egw in egw_list:
            print("== Step 3:    GW ENCUTGW:%7.1f  ==" % egw)
            step_3_gw(exact_dir,gw_dir+'_egw_'+str(int(egw)),opts.tag_xc,opts.encut,nbands=nbands,nedos=nedos,encutgw=egw,nomega=opts.nomega,\
                  lwannier=opts.lwannier,vasp_cmd=vasp_cmd,gw_mode=opts.gw_mode,f_metal=opts.f_metal,smear=[0,0.05])

        gap_list = []
        for egw in egw_list:
            os.chdir(gw_dir+'_egw_'+str(int(egw)))
            gap_list.append(vasp_anal_get_outcar('gap'))
            os.chdir('..')
        os.chdir('..')

        #write to output file
        if os.path.exists('Gap_ENCUTGW.dat'):
            copy2('Gap_ENCUTGW.dat','Gap_ENCUTGW.dat_old')
        with open('Gap_ENCUTGW.dat','w') as f:
            f.write("#ENCUT=%6s MNPW=%6s NBANDS=%6s\n" % (encut,mnpw,nbands))
            f.write("#ENCUTGW gap\n")
            for i in xrange(len(egw_list)):
                f.write("%5.1f %12.6f\n" % (egw_list[i],gap_list[i]))

# =====================================================

if __name__ == "__main__":
    Main(sys.argv)
