#!/usr/bin/env python
# coding=utf-8

# ====================================================
#     File Name : pc_constants.py
# Creation Date : 26-04-2018
#    Created By : Min-Ye Zhang
#       Contact : stevezhang@pku.edu.cn
# ====================================================

import scipy.constants as scc

# Pi = 3.14159265....
Pi = scc.pi

# ==================================================

# fine stucture constant
fsa = scc.alpha
# Avogadro constant, in mol^-1
NAv = scc.N_A
# Planck constant, in J s
h = scc.h
hbar = scc.hbar
# Planck constant, in eV s
h_ev = scc.physical_constants['Planck constant in eV s'][0]
hbar_ev = scc.physical_constants['Planck constant over 2 pi in eV s'][0]
# mass of electron, in kg
me = scc.m_e
# elementary charge, in C
e = scc.e
# absolute g factor of electron
gf = scc.physical_constants['electron g factor'][0] * -1.0
# speed of light in vacuum, in meter
cLS = scc.c
# Boltzman constant, in J K^-1
kb = scc.k

# ==================================================

# Bohr radius, in meter
abohr = scc.physical_constants['Bohr radius'][0]

# ==================================================
# Energy conversion
eV2kcal      = 23.06
Ha2eV        = scc.physical_constants['Hartree energy in eV'][0]
eV2Ha        = 1.0/Ha2eV
Ry2eV        = Ha2eV/2.0
eV2Ry        = 1.0/Ry2eV
eV2J         = e
eV2kJ        = eV2J / 1000
eV2K         = e/kb
# kJ/mol to eV/f.u.
kJpMol2eVpFU = 1000 /(eV2J * NAv)
