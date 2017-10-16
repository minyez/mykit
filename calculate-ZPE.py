from sys import argv

# calculate ZPE from DFPT calculation

for outcar in argv[1:]:
    with open(outcar) as f:
        lines = f.readlines()

    print "================================"

    n_modes = 0
    n_imagi = 0
    zpe = 0.0

    for x in lines:
        line=x.split()
#        print len(line)
        if len(line) >= 7:
#            print line
            if line[4] == "THz":
#               print "    frequency %s: %s" % (line[0],line[9])
               n_modes = n_modes + 1
               zpe = zpe + float(line[9])
            elif line[3] == "THz":
#               print "(i) frequency %s: %s" % (line[0],line[8])
               n_modes = n_modes + 1
               n_imagi = n_imagi + 1
    zpe = zpe/2.0
    print " Zero-point energy: %s \n\
          # modes : %i \n\
  # imaginary freq: %i" % (zpe,n_modes,n_imagi)

    print "================================"


