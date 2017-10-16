#!/bin/env python
# get k-point-averaged band gap
# need files: IBZKPT, OUTCAR


klist = []  # [kx,ky,kz,weight]
gap = [] # [Eg]
VB = []
CB = []


with open("IBZKPT","r") as f:
    kptlines = f.readlines()
    nk = int(kptlines[1])
    for x in kptlines[3:3+nk]:
        klist.append(x.split())

with open("OUTCAR","r") as f:
    outlines = f.readlines()

nelec = 0
n = -1

for line in outlines:
    x = line.split()
    n += 1
    if len(x) > 0:
        if x[0] == "NELECT": # len(x) > 0 and 
             nelec = int(float(x[2]))
             nband = nelec/2
        if x[0] == "k-point" and len(x) == 6:
            HOMO,hocc = float(outlines[n+nband+1].split()[1]), outlines[n+nband+1].split()[2]
            LUMO,locc = float(outlines[n+nband+2].split()[1]), outlines[n+nband+2].split()[2]
            VB.append(HOMO)
            CB.append(LUMO)
#            print x[1],HOMO,hocc,LUMO,locc

kavgap = 0.0
VBM = VB[0]
CBM = CB[0]
k_VBM = [float(x) for x in klist[0][0:3]]
k_CBM = [float(x) for x in klist[0][0:3]]
kpts_tot = 0

for i in xrange(len(klist)):
    kpts = int(klist[i][3])
#    print i+1,kpts
    kpts_tot = kpts_tot + kpts
    kavgap = kavgap + ( CB[i] - VB[i] ) * kpts
    if VB[i] > VBM:
        VBM = VB[i]
        k_VBM = [float(x) for x in klist[i][0:3]]
    if CB[i] < CBM:
        CBM = CB[i]
        k_CBM = [float(x) for x in klist[i][0:3]]

print "CBM = %8.4f at (%6.4f, %6.4f, %6.4f)" % (CBM, k_CBM[0], k_CBM[1], k_CBM[2])
print "VBM = %8.4f at (%6.4f, %6.4f, %6.4f)" % (VBM, k_VBM[0], k_VBM[1], k_VBM[2])
print "Eg = %f" % (CBM-VBM)
print "K-averaged band gap = %f" % (kavgap/kpts_tot)
