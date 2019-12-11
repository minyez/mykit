#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from argparse import ArgumentParser
from mykit.wien2k.outputs import Scf


def pw_anal_scf():
    """Analyse the SCF file, e.g. energies, band gap, etc

    Currently only print out energy in each iteratio to 'Ener_Ite.dat' file.
    """
    parser = ArgumentParser(pw_anal_scf.__doc__)
    parser.add_argument("-f", dest="scf_file", default=None, help="name of SCF file")
    args = parser.parse_args()

    if args.scf_file is not None:
        scf = Scf(args.scf_file)
    else:
        scf = Scf()

    enes = []
    vols = []
    ites = []
    for ite, vals in scf.scf.items():
        ites.append(ite)
        enes.append(vals["ENE"])
        vols.append(vals["VOL"])
    # write energy of each iteration
    with open("Ener_Ite.dat", 'w') as h:
        for i, ite in enumerate(ites):
            print(ite, enes[i], file=h)


if __name__ == "__main__":
    pw_anal_scf()

