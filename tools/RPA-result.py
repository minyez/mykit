#!/usr/bin/python

from sys import argv

for x in argv[1:]:
    n = 0
    v = []
    eexx = []
    ecor = []
    etot = []
    with open(x) as f:
        lines = f.readlines()
    for i in lines:
        data=i.split()
#        print data
        n = n + 1
        v.append(data[1])
        eexx.append(data[3])
        ecor.append(data[4])
        etot.append(data[2])
    for i in [ v,eexx,ecor,etot]: print '    '.join(i)
