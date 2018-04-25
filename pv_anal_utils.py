#!/usr/bin/env python
#
# This script provide basic analysis tools for vasp output
#
from __future__ import print_function
import sys, os, shutil, copy
import subprocess as sp
import numpy as np
from pv_classes import vasp_read_xml
from scipy.optimize import curve_fit

# ====================================================
def vasp_anal_get_outcar(keyin,index=-1,outcar='OUTCAR'):
    '''
    return the value of key in outcar.
    '''
    key = keyin.lower()
# maximum number of plane wave, i.e. the maximum size of representation matrix
    if key=='mnpw':
        mnpw = sp.check_output("awk '/maximum number of/ {print $5}' %s | tail -1" % outcar,shell=True)
        return int(mnpw)
# NBANDS
    if key in ['nb', 'nbands']:
        nb = sp.check_output("awk '/NBANDS/ {print $15}' %s | head -1" % outcar,shell=True)
        return nb
# ENCUT
    if key in ['encut']:
        encut = sp.check_output("awk '/ENCUT/ {print $3}' %s | head -1" % outcar,shell=True)
        encut = int(float(encut))
        return encut
# converged G.S. energy
    if key in ['ene', 'energy']:
        ene = sp.check_output("awk '/without/ {print $7}' %s | tail -1" % outcar,shell=True)
        ene = float(ene)
        return ene
# total number of irreducible k-points
    if key in ['nkp']:
        nkp = sp.check_output("awk '/NKPTS/ {print $4}' %s | head -1" % outcar,shell=True)
        nkp = int(nkp)
        return nkp
# band gap. Use vaspxml class to obtain the value
    if key in ['gap','eg']:
        vaspxml = vasp_read_xml()
        gap = vaspxml.get_gap()
        return gap
    if key in ['efermi', 'e-fermi', 'fermi']:
        efermi = float(sp.check_output("awk '/E-fermi/ {print $3}' %s | tail -1" % outcar, shell=True))
        return efermi


def vasp_anal_get_enmax(potcar='POTCAR'):
    '''
    Get the largest ENMAX in the potcar file (default POTCAR)
    '''
    if not os.path.exists(potcar):
        print(" Error: POTCAR file not found: %s" % potcar)
        print(" Exit.")
        sys.exit(1)

    # sort in the increasing order, get the last value
    enmax_str = sp.check_output("awk '/ENMAX/ {print $3}' %s | sort | tail -1" % potcar,shell=True)
    # for example, ENMAX = 400.000;
    enmax = float(enmax_str.split(';')[0])
    return enmax


def vasp_anal_get_Ener_Vol(datatype='DFT'):
    '''
    Get calculated total energy and volume
    '''

# ====================================================

def vasp_anal_fit_EOS(name_ifile='Ener_Vol',nfu=1,eostype='BM',fixBp=False):
    '''
    Fit energy-volume data (Ener_Vol) with equation of states
    and give the fitting parameters
    Still need test

    :RETURN: fitting parameters, E in eV, V in A3, B in GPa
    BM: E0, V0, B0, B'
    '''
    def BMEOS(V,E0,V0,B0,Bp):
        return E0+9.0/16.0*V0*B0*((np.power(V0/V,2.0/3.0)-1.0)**3*Bp+(np.power(V0/V,2.0/3.0)-1.0)**2*(6.0-4.0*np.power(V0/V,2.0/3.0)))

# read Ener_Vol. n is the number of formula unit, default 1
    def read_Ener_Vol(filename,n):
        with open(filename) as f:
            lines = f.readlines()
        i = 0
        vol = [ ]
        ene = [ ]
        while i < len(lines):
            flag = 1
            line = lines[i].split()
    # skip commented line
            if line[0].startswith("#"):
                i +=1
                continue
    # avoid duplicates
            for v in vol:
                if v == float(line[1].strip()):
                    flag = 0
                    break
            if flag:
                vol.append(float(line[1].strip())/float(n))
                ene.append(float(line[2].strip())/float(n))
            i += 1
        data = [ np.array(vol), np.array(ene) ]
        return data

    data = read_Ener_Vol(name_ifile,nfu)

    if eostype == 'BM':
    # initialize fitting parameter
        E0 = min(data[1])
        bottom = np.where(data[1] == E0)
        V0 = data[0][bottom][0]
        B0 = 1.5
        if isinstance(fixBp,bool):
            if fixBp is True:
                print('fixBp not specified. Will not fix Bp')
            Bp = 4.0
            popt, pcov = curve_fit(BMEOS, data[0],data[1], p0 = [E0,V0,B0,Bp])
            Bp = popt[3]
        else:
            try:
                Bp = float(fitBp)
            except ValueError:
                raise "B' should be fixed with a float value"
            Bp = args.f
            popt, pcov = curve_fit(lambda V,E0,V0,B0: BMEOS(V,E0,V0,B0,Bp), data[0],data[1], p0 = [E0,V0,B0])
        opt_para = [popt[0], popt[1], popt[2], popt[3]]
# need to modify to output the errors
# return a list containing optimized parameters
    return opt_para

# ====================================================

def vasp_anal_read_eigen(spinpolarzied=False,debug=False):
    '''
    Read the EIGENVAL file and return the band_structure file

    # Currently it only supports reading EIGENVAL file
    # Return a list containing nkp members, each member being a list in a form of [[kx,ky,kz],energies of nband]

    '''
    f = open('EIGENVAL','r')

    lines = f.readlines()
    f.close()

    band_info = [int(x) for x in lines[5].split()]
    nelec, kpts_tot, nbands = band_info[0],band_info[1],band_info[2]

    band_struct = []
    band_struct.append(band_info)

    line_st = 7

    iline = line_st
    while iline < len(lines):
        kpt_band = []
        # add k-point information
        kpt_band.append(lines[iline].split())
        # add band energy
        for i in xrange(1,nbands+1):
            kpt_band.append(float(lines[iline+i].split()[1]))
        band_struct.append(kpt_band)
        iline += nbands + 2

    return band_struct

# ====================================================

def vasp_anal_get_BM_info(debug=False):
    '''
    Get irreducible k-points and HOMO-LUMO at each kpoint

    Return: k-points, top VB, bottom CB
    Need files: , or IBZKPT, OUTCAR
    '''

    klist = []  # [kx,ky,kz,weight]
    VB = []
    CB = []

    try:
        with open("EIGENVAL",'r') as f:
            lines = f.readlines()
    except:
        print("EIGENVAL file does not exist. Use IBZKPT and OUTCAR instead.")
    else:
        print("Read from EIGENVAL file.")
        band_info = [int(x) for x in lines[5].split()]
        nelec, kpts_tot, nbands = band_info[0],band_info[1],band_info[2]

        klist = []
        VB = []
        CB = []
        for k in xrange(kpts_tot):
    # get the full kpoint list: kx,ky,kz,weight(normalized)
            klist.append(lines[7+k*(nbands+2)].split())
    # get the valence band extreme
            VB.append(float(lines[7+k*(nbands+2)+nelec/2].split()[1]))
    # get the conduction band extreme
            CB.append(float(lines[7+k*(nbands+2)+nelec/2+1].split()[1]))
        return [klist, VB, CB]

    with open("IBZKPT","r") as f:
        kptlines = f.readlines()
        nk = int(kptlines[1])
        for x in kptlines[3:3+nk]:
#            kpoint = [float(y) for y in x.split()]
           # the weight, i.e. x.split()[3], is nor normalized
            klist.append(x.split())

    with open("OUTCAR","r") as f:
        outlines = f.readlines()

    nelec = 0

# find the last electronic iteration
    n_last_elec = int(sp.check_output("grep -n Iteration OUTCAR | tail -1",shell=True).split(":")[0]) - 2
    if debug: print(n_last_elec)

    n = 0
    while n < len(outlines):
        x = outlines[n].split()
        if len(x) > 0:
            if x[0] == "NELECT": # len(x) > 0 and
                 nelec = int(float(x[2]))
                 nband = nelec/2
                 # search the number of electron, then jump to the last iteration
                 n = n_last_elec - 1
                 if debug: print("Number of electrons: %i" % nelec)
            if x[0] == "k-point" and len(x) == 6:
            # occ is occupation, check for debug
                HOMO,hocc = float(outlines[n+nband+1].split()[1]), outlines[n+nband+1].split()[2]
                LUMO,locc = float(outlines[n+nband+2].split()[1]), outlines[n+nband+2].split()[2]
                VB.append(HOMO)
                CB.append(LUMO)
                if debug:
                    print(n,x[1],HOMO,hocc,LUMO,locc,x[3:6])
        n += 1
# return band extreme info
    return [ klist, VB, CB ]

# ====================================================

def vasp_anal_get_fund_gap(band_extreme_info,debug=False):   # deprecated
    '''
    Get fundamental band gap
    '''
    klist,VB,CB = band_extreme_info[0],\
                  band_extreme_info[1],\
                  band_extreme_info[2]

    VBM = max(VB)
    VBM_k_index = VB.index(VBM)
    VBM_k = [ float(x) for x in klist[VBM_k_index][0:3]]
    CBM = min(CB)
    CBM_k_index = CB.index(CBM)
    CBM_k = [ float(x) for x in klist[CBM_k_index][0:3]]
    kpts_tot = 0
    flag_metal = False
    if debug: print(len(klist))

    E_fd_gap = CBM - VBM
    print(" CBM = %8.4f eV   at (%7.4f,%7.4f,%7.4f)" % (CBM,CBM_k[0],CBM_k[1],CBM_k[2]))
    print(" VBM = %8.4f eV   at (%7.4f,%7.4f,%7.4f)" % (VBM,VBM_k[0],VBM_k[1],VBM_k[2]))
    print(" Fundamental gap Eg = %8.4f" % E_fd_gap)
    if E_fd_gap <= 0.0:
        flag_metal = True
    return E_fd_gap,flag_metal

# ====================================================

def vasp_anal_get_gap(band_struct,vb,cb,debug=False):
    '''
    Get minimal gap between two particular band
    '''

    nelec    = band_struct[0][0]
    nkp      = band_struct[0][1]
    band_max = band_struct[0][2]
    if debug:
        print(band_struct[0])

    VBM = max([band_struct[1+ikp][vb] for ikp in xrange(nkp)])
    VBM_k_index = [band_struct[1+ikp][vb] for ikp in xrange(nkp)].index(VBM)
    VBM_k = [ float(x) for x in band_struct[1+VBM_k_index][0][0:3]]
    E_gap_at_VBM = band_struct[1+VBM_k_index][cb] - band_struct[1+VBM_k_index][vb]
    if debug:
        print(VBM,VBM_k_index,VBM_k)

    CBM = min([band_struct[1+ikp][cb] for ikp in xrange(nkp)])
    CBM_k_index = [band_struct[1+ikp][cb] for ikp in xrange(nkp)].index(CBM)
    CBM_k = [ float(x) for x in band_struct[1+CBM_k_index][0][0:3]]
    E_gap_at_CBM = band_struct[1+CBM_k_index][cb] - band_struct[1+CBM_k_index][vb]
    if debug:
        print(CBM,CBM_k_index,CBM_k)

    E_gap = CBM - VBM
    print(" Band %3i:  BandMin = %8.4f eV   at (%7.4f,%7.4f,%7.4f)" % (cb, CBM, CBM_k[0], CBM_k[1], CBM_k[2]))
    print(" Band %3i:  BandMax = %8.4f eV   at (%7.4f,%7.4f,%7.4f)" % (vb, VBM, VBM_k[0], VBM_k[1], VBM_k[2]))
    print(" Eg(min) = %8.4f" % E_gap)
    if CBM_k_index != VBM_k_index:
        print(" Eg(VBM) = %8.4f" % E_gap_at_VBM)
        print(" Eg(CBM) = %8.4f" % E_gap_at_CBM)

# ====================================================

def vasp_anal_get_kavgap(band_struct,vb,cb,fix_k=-1,inv=False,debug=False):
    '''
    Get k-point-averaged band gap. Use fix_k to fix the k-point of a band, 0/1 for valence/conduction band
    Set inv to TRUE for the value of k-averaged inverse of band gap. This is for the analysis of polarization
    '''

    kpts_wtot = 0.0E0
    kavgap   = 0.0E0
    nelec    = band_struct[0][0]
    nkp      = band_struct[0][1]
    band_max = band_struct[0][2]

    VBM = max([band_struct[1+ikp][vb] for ikp in xrange(nkp)])
    VBM_k_index = [band_struct[1+ikp][vb] for ikp in xrange(nkp)].index(VBM)
    VBM_k = [ float(x) for x in band_struct[1+VBM_k_index][0][0:3]]

    CBM = min([band_struct[1+ikp][cb] for ikp in xrange(nkp)])
    CBM_k_index = [band_struct[1+ikp][cb] for ikp in xrange(nkp)].index(CBM)
    CBM_k = [ float(x) for x in band_struct[1+CBM_k_index][0][0:3]]

    kavgap = 0
    inv_kavgap = 0

    kpts_weigh = np.array([float(band_struct[i+1][0][3]) for i in xrange(nkp)])
    band_up = np.array([band_struct[1+i][cb] for i in xrange(nkp)])
    band_low = np.array([band_struct[1+i][vb] for i in xrange(nkp)])
    if fix_k == 0:
        band_low = np.fill(nkp,VBM)
    elif fix_k == 1:
        band_up = np.fill(nkp,CBM)
    else:
        pass

#   Sum of kavgap and inv_kavgap
    kavgap = np.inner( band_up - band_low, kpts_weigh )
    inv_kavgap = inv_kavgap + np.inner(np.reciprocal(band_up - band_low),kpts_weigh)

#   Average over kpoint
    kavgap = kavgap / np.sum(kpts_weigh)
    inv_kavgap = inv_kavgap / np.sum(kpts_weigh)
    if fix_k == -1:
        print("Mode -1: average direct band")
    elif fix_k == 0:
        print("Mode 0: average indirect band with transition to fixed VBM@(%7.4f,%7.4f,%7.4f)" % (VBM_k[0], VBM_k[1], VBM_k[2]))
    elif fix_k == 1:
        print("Mode 1: average indirect band with transition to fixed CBM@(%7.4f,%7.4f,%7.4f)" % (CBM_k[0], CBM_k[1], CBM_k[2]))
    else:
        pass
    if not inv:
        print(" E_{kav,g} (between band %3i and band %3i) = %8.4f eV" % (vb,cb,kavgap))
        return kavgap
    else:
        print(" E^{-1}_{kav,g} (between band %3i and band %3i) = %8.4f eV^{-1}" % (vb,cb,inv_kavgap))
        print(" Inverse of average inverse is %8.4f eV" % np.reciprocal(inv_kavgap))
        return inv_kavgap

# ====================================================
