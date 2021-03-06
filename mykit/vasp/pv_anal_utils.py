#!/usr/bin/env python
#
# This script provide basic analysis tools for vasp output
#

# from scipy.optimize import curve_fit



# ====================================================

# ====================================================

# def vasp_anal_read_eigen(spinpolarzied=False,debug=False):
#     '''
#     Read the EIGENVAL file and return the band_structure file

#     # Currently it only supports reading EIGENVAL file
#     # Return a list containing nkp members, each member being a list in a form of [[kx,ky,kz],energies of nband]

#     '''
#     f = open('EIGENVAL','r')

#     lines = f.readlines()
#     f.close()

#     band_info = [int(x) for x in lines[5].split()]
#     nelec, kpts_tot, nbands = band_info[0],band_info[1],band_info[2]

#     band_struct = []
#     band_struct.append(band_info)

#     line_st = 7

#     iline = line_st
#     while iline < len(lines):
#         kpt_band = []
#         # add k-point information
#         kpt_band.append(lines[iline].split())
#         # add band energy
#         for i in range(1,nbands+1):
#             kpt_band.append(float(lines[iline+i].split()[1]))
#         band_struct.append(kpt_band)
#         iline += nbands + 2

#     return band_struct

# ====================================================

# def vasp_anal_get_BM_info(debug=False):
#     '''
#     Get irreducible k-points and HOMO-LUMO at each kpoint

#     Return: k-points, top VB, bottom CB
#     Need files: , or IBZKPT, OUTCAR
#     '''

#     klist = []  # [kx,ky,kz,weight]
#     VB = []
#     CB = []

#     try:
#         with open("EIGENVAL",'r') as f:
#             lines = f.readlines()
#     except:
#         print("EIGENVAL file does not exist. Use IBZKPT and OUTCAR instead.")
#     else:
#         print("Read from EIGENVAL file.")
#         band_info = [int(x) for x in lines[5].split()]
#         nelec, kpts_tot, nbands = band_info[0],band_info[1],band_info[2]

#         klist = []
#         VB = []
#         CB = []
#         for k in range(kpts_tot):
#     # get the full kpoint list: kx,ky,kz,weight(normalized)
#             klist.append(lines[7+k*(nbands+2)].split())
#     # get the valence band extreme
#             VB.append(float(lines[7+k*(nbands+2)+nelec/2].split()[1]))
#     # get the conduction band extreme
#             CB.append(float(lines[7+k*(nbands+2)+nelec/2+1].split()[1]))
#         return [klist, VB, CB]

#     with open("IBZKPT","r") as f:
#         kptlines = f.readlines()
#         nk = int(kptlines[1])
#         for x in kptlines[3:3+nk]:
# #            kpoint = [float(y) for y in x.split()]
#            # the weight, i.e. x.split()[3], is nor normalized
#             klist.append(x.split())

#     with open("OUTCAR","r") as f:
#         outlines = f.readlines()

#     nelec = 0

# # find the last electronic iteration
#     n_last_elec = int(sp.check_output("grep -n Iteration OUTCAR | tail -1",shell=True).split(":")[0]) - 2
#     if debug: print(n_last_elec)

#     n = 0
#     while n < len(outlines):
#         x = outlines[n].split()
#         if len(x) > 0:
#             if x[0] == "NELECT": # len(x) > 0 and
#                  nelec = int(float(x[2]))
#                  nband = nelec/2
#                  # search the number of electron, then jump to the last iteration
#                  n = n_last_elec - 1
#                  if debug: print("Number of electrons: %i" % nelec)
#             if x[0] == "k-point" and len(x) == 6:
#             # occ is occupation, check for debug
#                 HOMO,hocc = float(outlines[n+nband+1].split()[1]), outlines[n+nband+1].split()[2]
#                 LUMO,locc = float(outlines[n+nband+2].split()[1]), outlines[n+nband+2].split()[2]
#                 VB.append(HOMO)
#                 CB.append(LUMO)
#                 if debug:
#                     print(n,x[1],HOMO,hocc,LUMO,locc,x[3:6])
#         n += 1
# # return band extreme info
#     return [ klist, VB, CB ]

# ====================================================

# def vasp_anal_get_gap(band_struct,vb,cb,debug=False):
#     '''
#     Get minimal gap between two particular band
#     '''

#     nelec    = band_struct[0][0]
#     nkp      = band_struct[0][1]
#     band_max = band_struct[0][2]
#     if debug:
#         print(band_struct[0])

#     VBM = max([band_struct[1+ikp][vb] for ikp in range(nkp)])
#     VBM_k_index = [band_struct[1+ikp][vb] for ikp in range(nkp)].index(VBM)
#     VBM_k = [ float(x) for x in band_struct[1+VBM_k_index][0][0:3]]
#     E_gap_at_VBM = band_struct[1+VBM_k_index][cb] - band_struct[1+VBM_k_index][vb]
#     if debug:
#         print(VBM,VBM_k_index,VBM_k)

#     CBM = min([band_struct[1+ikp][cb] for ikp in range(nkp)])
#     CBM_k_index = [band_struct[1+ikp][cb] for ikp in range(nkp)].index(CBM)
#     CBM_k = [ float(x) for x in band_struct[1+CBM_k_index][0][0:3]]
#     E_gap_at_CBM = band_struct[1+CBM_k_index][cb] - band_struct[1+CBM_k_index][vb]
#     if debug:
#         print(CBM,CBM_k_index,CBM_k)

#     E_gap = CBM - VBM
#     print(" Band %3i:  BandMin = %8.4f eV   at (%7.4f,%7.4f,%7.4f)" % (cb, CBM, CBM_k[0], CBM_k[1], CBM_k[2]))
#     print(" Band %3i:  BandMax = %8.4f eV   at (%7.4f,%7.4f,%7.4f)" % (vb, VBM, VBM_k[0], VBM_k[1], VBM_k[2]))
#     print(" Eg(min) = %8.4f" % E_gap)
#     if CBM_k_index != VBM_k_index:
#         print(" Eg(VBM) = %8.4f" % E_gap_at_VBM)
#         print(" Eg(CBM) = %8.4f" % E_gap_at_CBM)

# ====================================================

# def vasp_anal_get_kavgap(band_struct,vb,cb,fix_k=-1,inv=False,debug=False):
#     '''
#     Get k-point-averaged band gap. Use fix_k to fix the k-point of a band, 0/1 for valence/conduction band
#     Set inv to TRUE for the value of k-averaged inverse of band gap. This is for the analysis of polarization
#     '''

#     kpts_wtot = 0.0E0
#     kavgap   = 0.0E0
#     nelec    = band_struct[0][0]
#     nkp      = band_struct[0][1]
#     band_max = band_struct[0][2]

#     VBM = max([band_struct[1+ikp][vb] for ikp in range(nkp)])
#     VBM_k_index = [band_struct[1+ikp][vb] for ikp in range(nkp)].index(VBM)
#     VBM_k = [ float(x) for x in band_struct[1+VBM_k_index][0][0:3]]

#     CBM = min([band_struct[1+ikp][cb] for ikp in range(nkp)])
#     CBM_k_index = [band_struct[1+ikp][cb] for ikp in range(nkp)].index(CBM)
#     CBM_k = [ float(x) for x in band_struct[1+CBM_k_index][0][0:3]]

#     kavgap = 0
#     inv_kavgap = 0

#     kpts_weigh = np.array([float(band_struct[i+1][0][3]) for i in range(nkp)])
#     band_up = np.array([band_struct[1+i][cb] for i in range(nkp)])
#     band_low = np.array([band_struct[1+i][vb] for i in range(nkp)])
#     if fix_k == 0:
#         band_low = np.fill(nkp,VBM)
#     elif fix_k == 1:
#         band_up = np.fill(nkp,CBM)
#     else:
#         pass

#   Sum of kavgap and inv_kavgap
#     kavgap = np.inner( band_up - band_low, kpts_weigh )
#     inv_kavgap = inv_kavgap + np.inner(np.reciprocal(band_up - band_low),kpts_weigh)

# #   Average over kpoint
#     kavgap = kavgap / np.sum(kpts_weigh)
#     inv_kavgap = inv_kavgap / np.sum(kpts_weigh)
#     if fix_k == -1:
#         print("Mode -1: average direct band")
#     elif fix_k == 0:
#         print("Mode 0: average indirect band with transition to fixed VBM@(%7.4f,%7.4f,%7.4f)" % (VBM_k[0], VBM_k[1], VBM_k[2]))
#     elif fix_k == 1:
#         print("Mode 1: average indirect band with transition to fixed CBM@(%7.4f,%7.4f,%7.4f)" % (CBM_k[0], CBM_k[1], CBM_k[2]))
#     else:
#         pass
#     if not inv:
#         print(" E_{kav,g} (between band %3i and band %3i) = %8.4f eV" % (vb,cb,kavgap))
#         return kavgap
#     else:
#         print(" E^{-1}_{kav,g} (between band %3i and band %3i) = %8.4f eV^{-1}" % (vb,cb,inv_kavgap))
#         print(" Inverse of average inverse is %8.4f eV" % np.reciprocal(inv_kavgap))
#         return inv_kavgap

# ====================================================
