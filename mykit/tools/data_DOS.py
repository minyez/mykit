#!/usr/bin/python
import math
pi = 3.14159265
sqrt_pi = (pi)**0.5
with open('mesh.yaml') as f:
    lines = f.readlines()
nqp = int(lines[1].split(":")[-1].strip())
nf = int(lines[2].split(":")[-1].strip()) * 3
i = 8
j = 0
n = 2 * nf + 4
sigma = 1E-1
weight = []
freq = []
prob = []
while (i+n-1) <= len(lines):
    x = [float(lines[i+1].split(":")[-1].strip())] * nf
    weight.extend(x)
    while j <= (nf-1):
        freq.append(float(lines[i+4+j*2].split(":")[-1].strip()))
	j += 1
    i += n
    j = 0
mx = max(freq)
mn = min(freq)
sta = mn - 3.0*sigma
end = mx + 3.0*sigma
pos = sta
Z = 0.0		# Z is the normalization factor
tempp = 0.0
while pos >= sta and pos <= end:
	while j <= (len(freq)-1):
	    if (freq[j] >= (pos-3.0*sigma) and freq[j] <=(pos+3.0*sigma)):
		tempp = tempp + weight[j]*math.exp(-(pos-freq[j])**2/sigma**2)/sqrt_pi/sigma
	    j += 1
	prob.append([pos,tempp])
	Z = Z + tempp
	tempp = 0.0
	j = 0
	pos = pos + sigma
dataout = open("dosdist.dat","w")
for x in prob
    x[0] = str(x[0])
    x[1] = str(x[1]/Z)
    dataout.write(" ".join(x)+'\n')
dataout.close()
