#!/usr/bin/env python
# coding=utf-8

# ====================================================
#
#     File Name :
# Creation Date : 2017-05-01
# Last Modified : Thu 30 Nov 2017 09:15:02 AM CST
#    Created By : Min-Ye Zhang
#       Contact : stevezhang@pku.edu.cn
#       Purpose : This script offers utilities for initializing and running
#                 tasks including standard DFT, LDA/GGA+U, hybrid functional
#                 , GW and RPA calculations as well as some basic io utilities
#                 Special thanks to Prof. Hong Jiang and Dr. Feng Wu
#
# ====================================================

import sys, re, os, shutil, copy
import subprocess as sp
import commands, string
from argparse import ArgumentParser
from math import sqrt
#from io_utils import *

# ====================================================
# need bugfix
def vasp_io_get_tag_value(tag,n_val=1,ifile='INCAR',eq='=',debug=False):
    '''
    Get the value correspondent to the tag in the ifile (default INCAR). The case of tag does not matter.
    :Return: an object whose type is correspondent to the tag
    '''
    tag_upcase = tag.upper()
    val = None
    if ifile == 'INCAR':
        with open(ifile,'r') as f:
            lines = f.readlines()
            for line in lines:
                x = [ y.strip() for y in re.split(r'[%s;]'%eq,line)]
                if debug: print x
                if '#' in x:
                    i_comment = x.index('#')
                    if debug: print i_comment
                    if not i_comment == 0:
                        del x[i_comment:]
                    else: # this is just a comment line
                        continue
                if tag in x: # usually tags in INCAR are uppercas
                    i_tag = x.index(tag)
                    if n_val == 1:
                        val = x[i_tag+1]
                    else:
                        val = x[i_tag+1:i_tag+n_val+1]

    return val

# ====================================================

def vasp_io_change_tag(name_ifile,tag,name_ofile='temp',new_val=None,n_val=1,replace=False,eq='=',backup=True):

# add, delete tag or change tag value in the name_ifile
# for now can change one parameter and single-value parameter. If new_val=None, the tag will be delete
# change MAGMOM LDA+U tags are not supported yet
    ifile=open(name_ifile,'r')
    ofile=open(name_ofile,'w')

    lines = ifile.readlines()
    flag_tag = 0
    length = len(lines)
    for i in xrange(length):
    # split the line with eq notation and semicolon
        line = [ x.strip() for x in re.split(r'[%s;]'%eq,lines[i])]
        if len(line) == 0:
            ofile.write(lines[i])
            continue
        elif len(line) >= 2 and line[0] != "#":
            if tag in line:
                flag_tag = 1
                j = line.index(tag)
                if str(new_val) == "None":
                    del line[j], line[j]
                else:
                    line[j+1] = new_val
                l = len(line)
                if l == 0: continue
                else:
                    taglist = ['']
                    k = 0
                    while k+1 < l:
                        if line[k] == "#" or line[k+1] == "#":
                            break
                        templist = [str(line[k]),eq,str(line[k+1]),";"]
                        taglist = taglist + templist
                        k = k + 2
                    newline = ' '.join(taglist)+'\n'
                    # the last ';' identifies the modified line
                    ofile.write(newline)
            else:
                ofile.write(lines[i])
        else:
            ofile.write(lines[i])

# if the tag is not found, add it to the bottom of file
    if flag_tag == 0 and str(new_val) != 'None':
        print '%s tag does not exist. Add it to %s' % (tag,name_ifile)
        ofile.write(' %s = %s\n'%(tag,new_val))

    ofile.close()
    ifile.close()

# if replace is set or name_ofile is unset or names of ifile and ofile are the same, then save the ifile as ifile_bak and rename the ofile as ifile
    if replace or name_ofile == 'temp' or name_ifile == name_ofile:
        if backup:
            os.rename(name_ifile,name_ifile+'_bak')
        os.rename(name_ofile,name_ifile)

# ====================================================

def vasp_io_set_XC_type(incar,tag_xc):
    '''
    Set the tags associate with approximate functional
    set GGA tag for (semi-)local functionals and other tags for LDA+U or hybrid functional calculations
    '''

    tag_PBE0 = ' GGA = PE\n ALGO = ALL\n LHFCALC = .TRUE.\n TIME = 0.4\n PRECFOCK = Normal\n'

    xc_dict = {
        'LDA':    ' GGA = CA\n',
        'PBE':    ' GGA = PE\n',
        'PBEsol': ' GGA = PS\n',
        'RPBE':   ' GGA = RP\n',
        'PBE0':   '%s' % tag_PBE0,
        'HSE06':  '%s HFSCREEN = 0.2\n' % tag_PBE0,
        'HF':  '%s AEXX = 1.0\n' % tag_PBE0,
        'SCAN': ' METAGGA = SCAN\n',
        }

    if tag_xc is not None:
        tag_xc = tag_xc.strip()

        tag = xc_dict.get(tag_xc,None)
        if not tag == None:
            incar.write("\n# XC tag: %s\n" % tag_xc)
            incar.write(tag)
    # if tag is none, XC will be set according to POTCAR

# ====================================================

def common_io_checkdir(dirname=None,create=True):
    '''
    check if dirname exists, or create it
    return: the full path of target directory
    '''
    dirname = dirname.strip()

    if (dirname is None or dirname.strip() == ""):
        dirname = os.getcwd()
    elif (not os.path.exists(dirname)) and create:
        os.mkdir(dirname.strip())
    return dirname

# ====================================================

def common_io_cleandir(dirname=None):
    '''
    check if dirname exists and is empty, or create it
    and make it contain no files
    return: the full path of target directory
    '''
    if (dirname is None or dirname.strip() == ""):
        dirname = os.getcwd()
    elif (not os.path.exists(dirname)):
        os.mkdir(dirname)
    elif (os.path.exists(dirname)):
        sp.call("rm -rf %s" % (dirname+"/*"),shell=True)
    return dirname

# ====================================================

def vasp_io_get_NPAR(np):

    npar = 1
    nsqrt_max = int(sqrt(np)) + 1
    nsqrt = nsqrt_max
    while nsqrt > 0:
        if nsqrt % 2 == 0 and np % nsqrt == 0:
            npar = nsqrt
            break
        nsqrt = nsqrt - 1
    return npar

# ====================================================

def vasp_write_wannier90():
    '''
    '''

# ====================================================

def vasp_write_kpoints_basic(nks,mode='G',sh=None,debug=False,f_slab=None):

    if os.path.isfile('KPOINTS'):
        os.rename('KPOINTS','KPOINTS_old')

# if nks is put as a string
    if debug: print nks
    try:
        string = nks + 'str'
        nks = [ x.strip() for x in re.split(r'[,\[\]]',nks)]
        nks = [ int(x) for x in ' '.join(nks).split()]

# nks is input as a list of integer
    except TypeError:
        pass

# if slab model is used. f_slab = 1|2|3
    if f_slab is None:
        pass
    else:
        if f_slab in [1,2,3]:
            nks[f_slab-1] = 1

    if debug: print nks
    if mode == 'G':
        ofile = open('KPOINTS','w')
        ofile.write("K-Points\n")
        ofile.write("0 \n")
        ofile.write("Gamma \n")
        ofile.write("%d %d %d\n"%(nks[0],nks[1],nks[2]))
        if sh is None:
          ofile.write("0 0 0\n")
        else:
          ofile.write("%8.4f %8.4f %8.4f\n"%(sh[0],sh[1],sh[2]))

    elif mode == 'M':
        ofile = open('KPOINTS','w')
        ofile.write("K-Points\n")
        ofile.write("0 \n")
        ofile.write("Monkhorst-Pack \n")
        ofile.write("%d %d %d\n"%(nks[0],nks[1],nks[2]))
        if sh is None:
          ofile.write("0 0 0\n")
        else:
          ofile.write("%8.4f %8.4f %8.4f\n"%(sh[0],sh[1],sh[2]))
    ofile.close()


# ====================================================

def vasp_write_inputs_DOS_calc(ifile="INCAR"):
    '''
    Prepare the INCAR file for DFT calculation
    '''
    in_dos = 'INCAR_DOS'

    vasp_io_change_tag(ifile,'ISMEAR',in_dos,new_val='-5')
    vasp_io_change_tag(in_dos,'SIGMA')
    vasp_io_change_tag(in_dos,'ISTART')
    vasp_io_change_tag(in_dos,'LORBIT',new_val='11')
    vasp_io_change_tag(in_dos,'ICHARG',new_val='11')

# ====================================================

def vasp_write_incar_minimal_elec(incar,tag_xc,\
                encut=0,ediff="1E-8",npar=1,\
                mode_smear=None,wfrestart=False,chrestart=2,spin=1):
    '''
    write the minimum INCAR file
    '''
    incar.write("# Minimal electronic part\n")
    if not wfrestart:
        incar.write(" ISTART = 0\n")
        if chrestart != 2:
            incar.write(" ICHARG = %s\n" % chrestart)
    else:
        incar.write(" ISTART = %s\n" % wfrestart )
    if not encut == 0:
        incar.write(" ENCUT = %d\n" % encut)
        incar.write(" PREC = Accurate\n")
    else:
        incar.write(" PREC = Accurate # ENMAX will be used\n")
    incar.write(" LREAL = .False.\n")
    incar.write(" EDIFF = %s\n" % ediff)
    incar.write(" ISPIN = %s\n" % spin)

# mode_smear indicates the smearing setting and write it in INCAR.
# Should be None, a list like [ISMEAR, SIGMA], or an integer
#   None: write no smearing info (ISMEAR, SIGMA)
#      0: semiconductor smearing info, i.e. (0, 0.05)
#   1, 2: metal smearing info, i.e. (1, 0.2)
#     -1: Fermi smearing info, i.e. (-1, 0.05)
#     -5: Bl\:ochl smearing info, i.e. (-5) for bandstructure

    if not mode_smear == None:
        try:
            mode_smear = int(mode_smear)
            if mode_smear in [1,2]: incar.write(" ISMEAR = %d ; SIGMA = 0.2\n" % mode_smear)
            if mode_smear == -5: incar.write(" ISMEAR = %d\n" % mode_smear)
            if mode_smear == 0: incar.write(" ISMEAR = 0 ; SIGMA = 0.05\n")
            if mode_smear == -1: incar.write(" ISMEAR = 1 ; SIGMA = 0.05\n")
        except TypeError: # when a list parsed
            if len(mode_smear) == 2:
                incar.write(" ISMEAR = %s ; SIGMA = %s\n" % (mode_smear[0],mode_smear[1]))
            else: # reject other setting and will use default
                print "Invalid smearing setting, pass."

    if npar != 1:
        incar.write(" NPAR = %s\n" % npar)

# set XC functional
    vasp_io_set_XC_type(incar,tag_xc)

# ====================================================

def vasp_write_incar_LDAU(incar,tag_xc,encut,ediff="1E-8",mode=None):
    '''
    Under implementation !!!
    Write INCAR for LDA/GGA with Hubbard-U correction
    '''

# ====================================================

def vasp_write_incar_exact(incar,tag_xc,encut,nb=100,ediff="1E-8",npar=1,\
                           mode_smear=None,loptics=True,lwannier=False):
    '''
    Write INCAR for exact diagonalization for a particular number of nband. This is a restart calculation, default from WAVECAR.
    '''
    vasp_write_incar_minimal_elec(incar,tag_xc,encut,ediff=ediff,npar=npar,mode_smear=mode_smear,wfrestart=1)
    incar.write("\n# Exact Diag.\n")
    incar.write(" NELM   = 1\n")
    incar.write(" ALGO   = Exact\n")
    incar.write(" NBANDS = %s\n" % nb)
# For RPA calculation of system with a gap, WAVEDER is required to include longwave contribution. Set LOPTICS to generate it.
    if loptics: incar.write(" LOPTICS = .TRUE.\n")
# For GW and HF bandstructure calculation with WANNIER90
    if lwannier:
        incar.write(" LWANNIER90 = .TRUE.\n")
        vasp_write_wannier90()

# ====================================================

def vasp_write_incar_HF_coarse(incar,type_hf="HSE06",ediff="1E-6",npar=1,mode_smear=None):
    '''
    Under implementation !!!
    Write INCAR for hybrid functional with reduced k-mesh and a coarse FFT grid in HF routines (PRECFOCK = Fase). Note that for an accurate calculation, hybrid functional with original k-mesh and at least PRECFOCK = Normal is necessary.
    Supported HF: HSE06, PBE0
    '''

# ====================================================

def vasp_vaspcmd_zmy(np=1,mpitype="mpirun",vasp_path="vasp"):
    """
    Setting vasp command on different platforms
    mpitype indicates the type of MPI executives, whether use mpirun or yhrun
    """
# define execution type of vasp
    if vasp_path == "vasp":
        try:
            which_vasp = str(sp.check_output("which vasp",stderr=sp.STDOUT,shell=True)).split('\n')[0]
            vasp_path = which_vasp
        except sp.CalledProcessError:
            try:
                which_vasp = str(sp.check_output("which vasp_std",stderr=sp.STDOUT,shell=True)).split('\n')[0]
                vasp_path = which_vasp
            except sp.CalledProcessError:
                print "Path for vasp not found. Need manual setting"
                vasp_path = "vasp"
                raise

    if np == 1:
        vasp_cmd = vasp_path
    else:
        if mpitype == "mpirun":
            vasp_cmd = "mpirun -np %d %s" % (np, vasp_path)
        elif mpitype == "yhrun":
            vasp_cmd = "yhrun -n %d %s" % (np, vasp_path)

    return vasp_path, vasp_cmd

# ====================================================

def vasp_vasprun_zmy(vasp_cmd,fout=None,ferr=None):
    """
    Run vasp. mpi_cmd and vasp_path are both strings
    """

    if fout is None:
        ofile = sp.PIPE
    else:
        ofile = open(fout,'w')

    if ferr is None:
        efile = sp.PIPE
    else:
        efile = open(ferr,'w')

    p=sp.Popen(vasp_cmd,stdout=ofile,stderr=efile,shell=True)
    p.wait()

    if not fout is None: ofile.close()
    if not ferr is None: efile.close()

# ====================================================

