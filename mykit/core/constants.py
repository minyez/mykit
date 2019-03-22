# coding=utf-8
'''Module that defines math and physical constants to use

Constants from scipy.constants, 2019-03-07
'''
#pylint: disable=bad-whitespace

PI = 3.141592653589793
# fine stucture constant
FSCA = 0.0072973525664
# Avogadro constant, in mol^-1
NAV = 6.022140857e+23
# Planck constant, in J s
H = 6.62607004e-34
HBAR = 1.0545718001391127e-34
# Planck constant, in eV s
HINEV = 4.135667662e-15
HBARINEV = 6.582119514e-16
# mass of electron, in kg
ME = 9.10938356e-31
# elementary charge, in C
E = 1.6021766208e-19
# absolute g factor of electron
GF = 2.00231930436182
# speed of light in vacuum, in meter
CLIGHT = 2.99792458e8
# Boltzman constant, in J K^-1
KB = 1.38064852e-23
# Bohr radius, in meter
AU = 5.2917721067e-11

# ==================================================
# Length conversion
AU2ANG = AU * 1.0E10
ANG2AU = 1.0 / AU2ANG
# Energy conversion
EV2KCAL = 23.060548
HA2EV = 27.21138602
EV2HA = 1.0 / HA2EV
RY2EV = HA2EV / 2.0
RY2HA = 0.5
EV2RY = 1.0 / RY2EV
EV2J = E
EV2KJ = EV2J / 1000.0
EV2K = E / KB
# kJ/mol to eV/f.u.
KJ_PER_MOL2EV_PER_FU = 1000 / (EV2J * NAV)
# eV/A^3 to GPa
EV_PER_ANG_CUB2GPA = E * 1e21
