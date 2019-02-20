# coding=utf-8
'''Module that defines math and physical constants to use
'''
#pylint: disable=bad-whitespace

import scipy.constants as scc

pi = scc.pi
# fine stucture constant
fsa = scc.alpha
# Avogadro constant, in mol^-1
nav = scc.N_A
# Planck constant, in J s
h = scc.h
hbar = scc.hbar
# Planck constant, in eV s
hInEv = scc.physical_constants['Planck constant in eV s'][0]
hbarInEv = scc.physical_constants['Planck constant over 2 pi in eV s'][0]
# mass of electron, in kg
me = scc.m_e
# elementary charge, in C
e = scc.e
# absolute g factor of electron
gf = scc.physical_constants['electron g factor'][0] * -1.0
# speed of light in vacuum, in meter
cLight = scc.c
# Boltzman constant, in J K^-1
kb = scc.k

# Bohr radius, in meter
au = scc.physical_constants['Bohr radius'][0]

# ==================================================
# Length conversion
au2ang  = au * 1.0E10
ang2au  = 1.0 / au2ang
# Energy conversion
ev2kcal = 23.060548
ha2ev   = scc.physical_constants['Hartree energy in eV'][0]
eV2ha   = 1.0 / ha2ev
ry2ev   = ha2ev / 2.0
ev2ry   = 1.0 / ry2ev
ev2j    = e
ev2kj   = ev2j / 1000.0
ev2k    = e / kb
# kJ/mol to eV/f.u.
kjPerMol2evPerFu = 1000 /(ev2j * nav)

