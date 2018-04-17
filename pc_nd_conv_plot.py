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
    #para_names = [x for x in df.columns]
    # Get the column names and the maximum value for each column
    # Here the fact that the calculation is more accurate with larger parameter is assumed.

    # For n_columns >= 5, it is not recommended
    if n_columns >= 5:
        raise ValueError(" data columns >= 5 will be crowded and NOT implemented YET. Remove some data.")
    if ytarget_column == 0:
        i_ytarget = n_columns - 1
    else:
        try:
            assert ytarget_column <= n_columns
            assert ytarget_column > 0
        except AssertionError:
            raise ValueError("-t value should be non-negative and less equal to -n")
        else:
            i_ytarget = ytarget_column - 1

    if xtarget_column == 0:
        i_xtarget = n_columns - 2
    else:
        try:
            assert xtarget_column <= n_columns
            assert xtarget_column > 0
        except AssertionError:
            raise ValueError("-t value should be non-negative and less equal to -n")
        else:
            i_xtarget = xtarget_column - 1

    para_names = []
    for i in range(n_columns):
        if i == i_xtarget or i == i_ytarget:
            continue
        para_names.append(df.columns[i])

    col_max = []
    for col in para_names:
        col_max.append(df[col].max())

    x_name = df.columns[i_xtarget]
    y_name = df.columns[i_ytarget]

    return n_columns, x_name, y_name, para_names, col_max

# ====================================================

def __init_fig_axs(n_columns, para_names, x_name, y_name):
  
    if n_columns == 3:
    # only 1 graph is sufficient for 2 convergence parameters, 
    # with one as the x-axis
        fig, axs = plt.subplots()
        axs[0].set_xlabel(x_name, size=12)
        axs[0].set_ylabel(y_name,size=12)

    elif n_columns == 4:
    # At least 2 graphs are required for 3 convergence parameters,
    # with one as the x-axis
        fig, axs = plt.subplots(1,2)
        axs[0].set_xlabel(x_name, size=12)
        axs[1].set_xlabel(x_name, size=12)
        axs[0].set_ylabel(y_name, size=12)

    return fig, axs

# ====================================================

def __init_fig_3d_axs(n_columns, para_names, x_name, y_name):
    pass

# ====================================================

def pc_nd_conv_plot(datafile, xtarget_column, ytarget_column, f_plot3d, figname):

    df_all = pd.read_table(datafile, delim_whitespace=True)

    n_columns, x_name, y_name, para_names, col_max = \
            __check_column_and_target(df_all, xtarget_column, ytarget_column)

    # Group the DataFrame by groupby method
    df_all_gpb = df_all.groupby(para_names)

    # TODO:
    #   if 3D plot is required, import necessary 3D plotting modules first
    if f_plot3d:
        fig, axs = __init_fig_3d_axs(n_columns, para_names, x_name, y_name)

    else:
        fig, axs = __init_fig_axs(n_columns, para_names, x_name, y_name)

        if n_columns == 3:
            for group in df_all_gpb.groups:
                gp_data = df_all_gpb.get_group(group)
                axs[0].plot(gp_data[x_name], gp_data[y_name], label="%s=%s" % (para_names[0], group[0]))

        if n_columns >= 4:
            for group in df_all_gpb.groups:
                for i in range(len(para_names)):
                # check the convergence of parameter para_names[i]
                # with the other parameters at the best, i.e. max
                    flag_best_other = True
                    for j in range(len(para_names)):
                        if j != i and group[j] != col_max[j]:
                            flag_best_other = False
                            break
                    if not flag_best_other:
                        continue

                    gp_data = df_all_gpb.get_group(group)
                    axs[i].plot(gp_data[x_name], gp_data[y_name], label="%s=%s" % (para_names[i], group[i]))

            # Generate the title string as the fixed parameters
            for i in range(len(para_names)):
                title_str_list = ['convergence w.r.t', para_names[i],'\n@ (']
                for j in range(len(para_names)):
                    if j != i:
                        title_str_list.append("%s = %s" % (para_names[j], col_max[j]))
                title_str_list.append(')')
                title_str = ' '.join(title_str_list)
                axs[i].set_title(title_str)

        for ax in axs:
            ax.legend(loc="upper left", shadow=True, fancybox=True)

    plt.show()

    if figname is not '':
        print("- Saving to %s" % figname)
        fig.savefig(figname, dpi=600)

    return

# ====================================================

def Main(ArgList):

    description = '''Visualize the data for an N-parameter convergence test. In general N is equal to 2 or 3.''' 
    parser = ArgumentParser(description=description)
    parser.add_argument(dest="datafile", metavar='file', type=str, nargs=1, help="The name of file storing the data. Better in CSV/Excel format and index is not necessary.")
    #parser.add_argument("-n", dest="N", type=int, default=0, help="the number of columns in the data file (except the index). Default 0 for automatic detection.")
    parser.add_argument("--yt", dest="ytarget_column", type=int, default=0, help="the index of column (>0) which contains the quantity to converge (y). Default is the last column.")
    parser.add_argument("--xt", dest="xtarget_column", type=int, default=0, help="the index of column (>0) which contains the direct test parameter (x). Default is the second to last column.")
    parser.add_argument("--plot3d", dest="f_plot3d", action="store_true", help="Flag to use 3D plots. NOT implemented YET.")
    parser.add_argument("--save", dest="figname", type=str, default='', help="File name (e.g. conv.png) to save the figure. The figure will not be saved unless this option is set other than ''.")

    # initialize options as 'opts'
    opts = parser.parse_args()
    datafile = opts.datafile[0]

    pc_nd_conv_plot(datafile, opts.xtarget_column, opts.ytarget_column, opts.f_plot3d, opts.figname)

# ==============================

if __name__ == "__main__":
    Main(sys.argv)

