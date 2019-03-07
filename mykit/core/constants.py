# coding=utf-8
'''Module that defines math and physical constants to use

Constants from scipy.constants, 2019-03-07
'''
#pylint: disable=bad-whitespace

pi = 3.141592653589793
# fine stucture constant
fsca = 0.0072973525664
# Avogadro constant, in mol^-1
nav = 6.022140857e+23
# Planck constant, in J s
h = 6.62607004e-34
hbar = 1.0545718001391127e-34
# Planck constant, in eV s
hInEv = 4.135667662e-15
hbarInEv = 6.582119514e-16
# mass of electron, in kg
me = 9.10938356e-31
# elementary charge, in C
e = 1.6021766208e-19
# absolute g factor of electron
gf = 2.00231930436182
# speed of light in vacuum, in meter
cLight = 2.99792458e8
# Boltzman constant, in J K^-1
kb = 1.38064852e-23
# Bohr radius, in meter
au = 5.2917721067e-11

# ==================================================
# Length conversion
au2ang  = au * 1.0E10
ang2au  = 1.0 / au2ang
# Energy conversion
ev2kcal = 23.060548
ha2ev   = 27.21138602
eV2ha   = 1.0 / ha2ev
ry2ev   = ha2ev / 2.0
ev2ry   = 1.0 / ry2ev
ev2j    = e
ev2kj   = ev2j / 1000.0
ev2k    = e / kb
# kJ/mol to eV/f.u.
kjPerMol2evPerFu = 1000 /(ev2j * nav)
# eV/A^3 to GPa
evPerAcub2gpa = e * 1e21
