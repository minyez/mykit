#!/usr/bin/env python
# coding=utf-8

# ====================================================
#     File Name : pc_gw_extrap.py
# Creation Date : 09-05-2018
#    Created By : Min-Ye Zhang
#       Contact : stevezhang@pku.edu.cn
# ====================================================

from __future__ import print_function
import sys
import numpy as np
from scipy.optimize import curve_fit
from argparse import ArgumentParser
import matplotlib.pyplot as plt

def __gw_gap_vs_nband(nband, num, denom, eg_inf):
    return num/(nband - denom) + eg_inf


def common_gw_gap_extrapolate(array_nband, array_gwgap):
    popt, pcov = curve_fit(__gw_gap_vs_nband, array_nband, array_gwgap, p0=[-10.0,-1.0,array_gwgap.max()])
    return popt


def __get_gap_data(filename):

    with open(filename,'r') as h_gb:
        datalines = h_gb.readlines()
    nband = []
    gwgap = []
    for line in datalines:
        if line.strip().startswith("#"):
            continue
        words = line.split()
        nband.append(words[0])
        gwgap.append(words[1])
    array_nband = np.array(nband, dtype='float64')
    array_gwgap = np.array(gwgap, dtype='float64')

    return array_nband, array_gwgap


def __plot_scatter(ax, array_nband, array_gwgap):

    ax.plot(array_nband, array_gwgap, 'ob')


def __plot_fit_curve(ax, params, nbandmin, nbandmax):

    num    = params[0]
    denom  = params[1]
    eg_inf = params[2]

    array_nband = np.linspace(nbandmin, nbandmax, 500)
    array_gwgap = __gw_gap_vs_nband(array_nband, num, denom, eg_inf)
    ax.axhline(eg_inf, c='green', linestyle='dashed', linewidth=2)
    ax.plot(array_nband, array_gwgap, '-g')


def Main(ArgList):

    description='''extrapolate the GW band gap to infinite NBANDS'''
    
    parser = ArgumentParser(description=description)
    
    parser.add_argument('filename', nargs='+', help="optional argument with undetermined number of paras")
    parser.add_argument("-p", dest="plot", action="store_true", help="flag for plot the data")
    parser.add_argument("-D", dest="debug", action="store_true", help="flag for debug mode")
    
    # initialize options as 'opts'
    opts = parser.parse_args()

    if opts.plot:
        fig, ax = plt.subplots(figsize=(8,8))

    for i in range(len(opts.filename)):
        filename = opts.filename[i]
        array_nband, array_gwgap = __get_gap_data(filename)
        params = common_gw_gap_extrapolate(array_nband, array_gwgap)
        print("File %2d: Converged GW gap: %8.4f" % (i+1, params[2]))
        if opts.debug:
            print(params)
        if opts.plot:
            ax.clear()
            __plot_scatter(ax, array_nband, array_gwgap)
            __plot_fit_curve(ax, params, array_nband.min()*0.95, array_nband.max()*1.05)
            plt.show()
            

# ==============================

if __name__ == "__main__":
    Main(sys.argv)

