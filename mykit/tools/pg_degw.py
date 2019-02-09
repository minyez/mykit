#!/usr/bin/env python
# coding=utf-8

# ====================================================
#     File Name : pg_degw.py
# Creation Date : 05-09-2018
#    Created By : Min-Ye Zhang
#       Contact : stevezhang@pku.edu.cn
# ====================================================

from __future__ import print_function, absolute_import
import sys
import os
from argparse import ArgumentParser
import numpy as np
import matplotlib.pyplot as plt

def __datafile_name(mainname, fnpre=None, fnsuf=None):
    '''
    Generate the filename with specific suffix and prefix of .dat file

    Parameters
    ----------
    mainname : str
    fnsuf : str
    fnpre : str

    Returns
    -------
    str : final name of the .dat file

    Examples
    --------
    >>>
    '''

    DataFile = mainname

    if fnpre is not None:
        DataFile = fnpre + '_' + DataFile
    if fnsuf is not None:
        DataFile = DataFile + '_' + fnsuf
    DataFile = DataFile + '.dat'

    return DataFile

# ==============================

def read_eqpeV(eqpeV_file, use_EKS_x, col=-4):
    '''
    Read eqpeV_file to get quasi-particle correction

    Parameters
    ----------
    eqpeV_file : str
        filename to extract the degw data
    use_EKS_x : bool
        flag for using Kohn-Sham energy level as x-axis
    col : int
        this param indicates the location of y data,
        default degw is in the 4 column to the last

    Returns
    -------
    int : number of valence bands
    list :
    list :

    Examples
    --------
    >>>
    '''


    degw_data = []
    x_data = []
    k_points = []
    vb = 0

    with open(eqpeV_file, 'r') as h_eqp:
        eqp_lines = h_eqp.readlines()

    bandl, bandh = eqp_lines[4].split()[-2:]
    bandl = int(bandl)
    bandh = int(bandh)
    nk = int(eqp_lines[3].split()[-1])
    nspin = int(eqp_lines[2].split()[-1])

    bandgw = bandh - bandl + 1
    
    for ik in range(nk):
        degw_data.append([])
        x_data.append([])
        kline = eqp_lines[10 + ik * (bandgw + 2)]
        kint = [int(kline[25 + i * 4:29 + i * 4]) for i in range(3)]
        kdiv = int(kline[44:48])
        k_points.append([float(x)/kdiv for x in kint])
        datalines = eqp_lines[12 + ik * (bandgw + 2):12 + ik * (bandgw + 2) + bandgw]
        for line in datalines:
            degw_data[-1].append(float(line.split()[col]))
            EKS = line.split()[2]
            iband = int(line.split()[1])
            # automatic find the index of top valence band
            # currently only for semiconductor and insulator,
            # since all kpoints have the same top valence band index
            if EKS == "0.000" and iband > vb:
                vb = iband
            if use_EKS_x:
                x_data[-1].append(float(EKS))
            else:
                x_data[-1].append(float(iband))

    # arrange the DEGWs of each band in corresponding list, i.e. transpose the data
    degw_data = np.transpose(np.array(degw_data))
    x_data = np.transpose(np.array(x_data))
    nvb = vb - bandl + 1

    return nvb, x_data, degw_data, k_points

# ==============================

def __export_data(x_data, degw_data, nvb, data_index, filename, fnpre, fnsuf):
    '''
    Export the degw data
    Parameters
    ----------
    x_data : ndarray
    degw_data : ndarray
    nvb : int
    data_index : int
    filename : str
    fnpre : str
    fnsuf : str

    Returns
    -------
    None

    Examples
    --------
    >>>
    '''
    # split data into valence and conduction part for coloring
    VB_data = [[], []]
    CB_data = [[], []]

    for i in range(len(x_data)):
        if i < nvb:
            data_region = VB_data
        else:
            data_region = CB_data
        for x, degw in zip(x_data[i], degw_data[i]):
            data_region[0].append(x)
            data_region[1].append(degw)

    bandtypes = {"VB": VB_data, "CB": CB_data}

    for btype in bandtypes.iterkeys():
        OutFile = __datafile_name("degw_%s_%02d" % (btype, data_index), fnpre, fnsuf)
        with open(OutFile, 'w') as h_out:
            h_out.write("#data from %s\n" % os.path.abspath(filename))
            h_out.write("#x-axis could be band index or KS energy level\n")
            h_out.write("#x degw(eV)\n")
            for x, degw in zip(bandtypes[btype][0], bandtypes[btype][1]):
                h_out.write("%7.3f  %9.3f\n" % (x, degw))

# ==============================

def __plot_degw(axs, use_EKS_x, f_compare, fnpre, fnsuf):
    '''
    Plot the degw with matplotlib. In maximum 2 sets of data will be plot

    Parameters
    ----------
    axs : matplotlib axis object
    use_EKS_x : bool
    f_compare : bool
    fnpre : str
    fnsuf : str

    Returns
    -------
    None

    Examples
    --------
    >>>
    '''
    if use_EKS_x:
        axs.set_xlabel("$\epsilon_{KS}$ (eV)", fontsize=14)
    else:
        axs.set_xlabel("Band index", fontsize=14)
    axs.set_ylabel("$\Delta\epsilon$ (eV)", fontsize=14)

    if f_compare:
        len_datafile = 2
    else:
        len_datafile = 1

    bandtypes = {"VB": 'r', "CB": 'b'}

    # check the range of x and y in data
    xmins = []
    xmaxs = []
    ymins = []
    ymaxs = []

    for i in range(len_datafile):
        data_index = i + 1
        for btype in bandtypes.iterkeys():
            x = []
            y = []
            DataFile = __datafile_name("degw_%s_%02d" % (btype, data_index), fnpre, fnsuf)
            with open(DataFile, 'r') as h_in:
                datalines = h_in.readlines()[3:]
                for line in datalines:
                    x.append(float(line.split()[0]))
                    y.append(float(line.split()[1]))
                xmins.append(min(x))
                xmaxs.append(max(x))
                ymins.append(min(y))
                ymaxs.append(max(y))
                if i:
                    axs.plot(x, y, 'o'+bandtypes[btype], markersize=10)
                else:
                    axs.plot(x, y, 'o'+bandtypes[btype], markersize=10, mfc="none")

        # for legend only
        if i:
            axs.plot([1000], [100], 'ok', markersize=10, label='data 2')
        else:
            axs.plot([1000], [100], 'ok', markersize=10, mfc="none", label='data 1')

    # plot zero
    ngrid = 200
    axs.set_xlim([min(xmins) - 2,  max(xmaxs) + 2])
    axs.set_ylim([min(ymins) - 0.5, max(ymaxs) + 0.5])
    zero_x = np.linspace(min(xmins ) -2, max(xmaxs) + 2, ngrid)
    zero_y = np.zeros(ngrid)
    axs.plot(zero_x, zero_y, linestyle="dashed", color="black")
    axs.legend()

# ==============================

def __Main(ArgList):

    description = '''Extract the self-energy correction ($\Delta\epsilon$,`degw`) \
            from GAP calculation (case.eqpeV_GW/GW0)'''

    parser = ArgumentParser(description=description)

    parser.add_argument(dest="filenames", nargs='+', \
            help="eqpeV file names. No more than 2 for comparison")
    parser.add_argument("-e", dest="f_EKS_x", action="store_true", \
            help="Flag to use KS energy level as x-axis")
    parser.add_argument("-p", dest="f_plot", action="store_true", \
            help="Flag to plot the first two data")
    parser.add_argument("--suf", dest="fnsuf", default=None, \
            help="Suffix to output files")
    parser.add_argument("--pre", dest="fnpre", default=None, \
            help="Prefix to output files")

    # initialize options as 'opts'
    opts = parser.parse_args()

    if opts.f_plot:
        fig, axs = plt.subplots(1, 1, figsize=(6, 6))

    for i in range(len(opts.filenames)):
        filename = opts.filenames[i]
        nvb, x_data, degw_data, k_points = read_eqpeV(filename, opts.f_EKS_x)
        __export_data(x_data, degw_data, nvb, i+1, filename, opts.fnpre, opts.fnsuf)
    # print VB and CD correction
        print(opts.filenames[i])
        print("%3s%-18s%12s%12s" % ("ik", " kvec  ", "DEGW_CB", "DEGWVB_VB"))
        for ik in range(len(k_points)):
            print("%3d%6.3f%6.3f%6.3f%12.6f%12.6f" % (ik+1, \
                    k_points[ik][0], k_points[ik][1], k_points[ik][2], \
                    degw_data[nvb][ik] , degw_data[nvb-1][ik]))

    if len(opts.filenames) > 1:
        f_compare = True
    else:
        f_compare = False
    if opts.f_plot:
        __plot_degw(axs, opts.f_EKS_x, f_compare, opts.fnpre, opts.fnsuf)
        plt.show()




# ==============================

if __name__ == "__main__":
    __Main(sys.argv)

