#!/usr/bin/env python

with open('thermal_properties.yaml') as f:
    lines = f.readlines()

i = 16
n = 6
lines2 = []
while (i+n-1) <= len(lines):
    lines2.append([x.split(":")[-1].strip() for x in lines[i:i+n-1]])
    i = i + n
dataout = open('thermdata.dat','w')
for line in lines2:
    dataout.write('  '.join(line)+'\n')
dataout.close()
