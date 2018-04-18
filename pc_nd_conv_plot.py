#!/usr/bin/env python
# coding=utf-8

# ====================================================
#
#     File Name : pc_nd_conv_plot.py
# Creation Date : 17-04-2018
#    Created By : Min-Ye Zhang
#       Contact : stevezhang@pku.edu.cn
#
# ====================================================

from __future__ import print_function
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser


def __check_column_and_target(df, xtarget_column, ytarget_column):

    n_columns = len(df.columns)
    
    # Get the column names and the maximum value for each column
    # Here the fact that the calculation is more accurate with larger parameter is assumed.

    # Not recommended to use for n_columns >= 7
    if n_columns >= 7:
        raise ValueError(" data columns >= 7 will be crowded and NOT implemented YET. Remove some data.")
    if ytarget_column == 0:
        i_ytarget = n_columns - 1
    else:
        try:
            assert ytarget_column <= n_columns
            assert ytarget_column > 0
        except AssertionError:
            raise ValueError("Invalid ytarget")
        else:
            i_ytarget = ytarget_column - 1

    if xtarget_column == 0:
        i_xtarget = n_columns - 2
    else:
        try:
            assert xtarget_column <= n_columns
            assert xtarget_column > 0
        except AssertionError:
            raise ValueError("Invalid xtarget")
        else:
            i_xtarget = xtarget_column - 1

    para_names = []
    for i in range(n_columns):
        if i == i_xtarget or i == i_ytarget:
            continue
        para_names.append(df.columns[i])

    para_max = []
    for col in para_names:
        para_max.append(df[col].max())

    x_name = df.columns[i_xtarget]
    y_name = df.columns[i_ytarget]

    return n_columns, x_name, y_name, para_names, para_max

# ====================================================

def __set_ax_linewidth(subplot_ax, linewidth=4):
    
    for axis in ['top','bottom','left','right']:
        subplot_ax.spines[axis].set_linewidth(linewidth)

    subplot_ax.tick_params(axis='both', which='major', length=linewidth*2, \
                           width=linewidth/2, direction='in')
    subplot_ax.tick_params(axis='both', which='minor', length=linewidth, \
                           width=linewidth/2, direction='in')

# ====================================================

def __init_fig_axs(n_columns, para_names, x_name, y_name):

    # N-1 graphs are required for N (n>=2) convergence parameters,
    # with the left one as the x-axis

    if n_columns == 3:
        fig, axs = plt.subplots(figsize=(8,8))
        axs.set_xlabel(x_name, size=12)
        axs.set_ylabel(y_name,size=12)
        __set_ax_linewidth(axs, 4)
    else:
        if n_columns == 4:
            fig, axs = plt.subplots(1,2, figsize=(12,8))
            axs[0].set_xlabel(x_name, size=12)
            axs[1].set_xlabel(x_name, size=12)
            axs[0].set_ylabel(y_name, size=12)
        if n_columns == 5:
            fig, axs = plt.subplots(1,3, figsize=(16,8))
            axs[0].set_xlabel(x_name, size=12)
            axs[1].set_xlabel(x_name, size=12)
            axs[2].set_xlabel(x_name, size=12)
            axs[0].set_ylabel(y_name, size=12)
        if n_columns == 6:
            fig, axs = plt.subplots(2,2, figsize=(12,12))
            #axs[:,:].set_xlabel(x_name, size=12)
            #axs[].set_xlabel(x_name, size=12)
            axs[0,0].set_ylabel(y_name, size=12)
            axs[1,0].set_ylabel(y_name, size=12)
            axs[1,0].set_xlabel(x_name, size=12)
            axs[1,1].set_xlabel(x_name, size=12)
        for ax in axs.flatten():
            __set_ax_linewidth(ax, 4)

    return fig, axs

# ====================================================

def __init_fig_3d_axs(n_columns, para_names, x_name, y_name):

    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure(figsize=(12,9))

    if n_columns == 3:
        axs = fig.add_subplot(111, projection='3d')
        axs.set_xlabel(para_names[0], size=12)
        axs.set_ylabel(x_name, size=12)
        axs.set_zlabel(y_name, size=12)

    return fig, axs

# ====================================================

def pc_nd_conv_plot(df_all, xtarget_column=0, ytarget_column=0, f_plot3d=False, \
                    figname='', preview=False, imgres=2):

    n_columns, x_name, y_name, para_names, para_max = \
            __check_column_and_target(df_all, xtarget_column, ytarget_column)

    # TODO:
    #   if 3D plot is required, import necessary 3D plotting modules first
    if f_plot3d:
        from matplotlib import cm

        fig, axs = __init_fig_3d_axs(n_columns, para_names, x_name, y_name)

        if n_columns == 3:
            p3d = axs.scatter(xs=df_all[para_names[0]], ys=df_all[x_name], zs=df_all[y_name], \
                        s=100, c=df_all[y_name], cmap=cm.coolwarm, marker='o', \
                        depthshade=False)
        else:
            raise ValueError("--plot3d has not been implemented for n_columns !=3. Sorry :(")

    else:
        # Group the DataFrame by groupby method
        df_all_gpb = df_all.groupby(para_names)

        fig, axs = __init_fig_axs(n_columns, para_names, x_name, y_name)

        if n_columns == 3:
            for group in df_all_gpb.groups:
                gp_data = df_all_gpb.get_group(group)
                axs.plot(gp_data[x_name], gp_data[y_name], 'o-', linewidth=2, \
                         label="%s=%s" % (para_names[0], group))
            axs.legend(loc="upper left", shadow=True, fancybox=True)

        if n_columns >= 4:
            for group in df_all_gpb.groups:
                for i in range(len(para_names)):
                # check the convergence of parameter para_names[i]
                # with the other parameters at the best, i.e. max
                    flag_best_other = True
                    for j in range(len(para_names)):
                        if j != i and group[j] != para_max[j]:
                            flag_best_other = False
                            break
                    if not flag_best_other:
                        continue

                    gp_data = df_all_gpb.get_group(group)
                    axs.flatten()[i].plot(gp_data[x_name], gp_data[y_name], 'o-', linewidth=2, \
                                label="%s=%s" % (para_names[i], group[i]))

            # Generate the title string as the fixed parameters
            for i in range(len(para_names)):
                title_str_list = ['convergence w.r.t', para_names[i],'\n@ (']
                for j in range(len(para_names)):
                    if j != i:
                        title_str_list.append("%s = %s" % (para_names[j], para_max[j]))
                title_str_list.append(')')
                title_str = ' '.join(title_str_list)
                axs.flatten()[i].set_title(title_str)

            for ax in axs.flatten():
                ax.legend(loc="upper left", shadow=True, fancybox=True)

    if preview:
        if f_plot3d:
            fig.colorbar(p3d)
        plt.show()

    if figname is not '':
        print("- Saving to %s" % figname)
        fig.savefig(figname, dpi=int(imgres)*150)

    return

# ====================================================

def Main(ArgList):

    description = '''Visualize the data for an N-parameter convergence test. In general N is equal to 2 or 3. Support up to 5.''' 
    parser = ArgumentParser(description=description)
    parser.add_argument(dest="datafile", metavar='file', type=str, nargs=1, help="The name of file storing the data. Better in CSV/Excel format and index is not necessary.")
    parser.add_argument("--xt", dest="xtarget_column", metavar="X", type=int, default=0, help="the index of column (>0) which contains the direct test parameter (x). Default is the second to last column.")
    parser.add_argument("--yt", dest="ytarget_column", metavar="Y", type=int, default=0, help="the index of column (>0) which contains the quantity to converge (y). Default is the last column.")
    parser.add_argument("--plot3d", dest="f_plot3d", action="store_true", help="Flag to use 3D plots. Support 2-parameter test only.")
    parser.add_argument("--save", dest="figname", type=str, default='', help="File name (e.g. conv.png) to save the figure. The figure will not be saved unless this option is set other than ''.")
    parser.add_argument("--res", dest="resolution", metavar='RES', type=int, default=2, help="Resolution of image, dpi = 150*RES. Default 2 (300 dpi).")

    # initialize options as 'opts'
    opts = parser.parse_args()
    datafile = opts.datafile[0]

    df_all = pd.read_table(datafile, delim_whitespace=True)

    pc_nd_conv_plot(df_all, opts.xtarget_column, opts.ytarget_column, opts.f_plot3d, opts.figname, \
                    True, opts.resolution)

# ==============================

if __name__ == "__main__":
    Main(sys.argv)

