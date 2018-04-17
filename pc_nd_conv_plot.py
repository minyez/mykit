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


def __check_column_and_target(n_columns, target_column):

    # For n_columns >= 5, it is not recommended
    if n_columns >= 5:
        raise ValueError(" data columns >= 5 will be crowded and NOT implemented YET. Remove some data.")
    if target_column == 0:
        i_target = n_columns - 1
    else:
        try:
            assert target_column <= n_columns
            assert target_column > 0
        except AssertionError:
            raise ValueError("-t value should be non-negative and less equal to -n")
        else:
            i_target = target_column - 1

    return i_target

# ====================================================

def __init_fig_axs(n_columns, col_names, t_col_name):
  
    if n_columns == 3:
    # only 1 graph is sufficient for 2 convergence parameters, 
    # with one as the x-axis
        fig, axs = plt.subplots()
        axs[0].set_xlabel(col_names[-1], size=12)
        axs[0].set_ylabel(t_col_name,size=12)

    elif n_columns == 4:
    # At least 2 graphs are required for 3 convergence parameters,
    # with one as the x-axis
        fig, axs = plt.subplots(1,2)
        axs[0].set_xlabel(col_names[-1], size=12)
        axs[1].set_xlabel(col_names[-1], size=12)
        axs[0].set_ylabel(t_col_name,size=12)

    return fig, axs

# ====================================================

def pc_nd_conv_plot(datafile, target_column, f_plot3d):

    df_all = pd.read_table(datafile, delim_whitespace=True)

    n_columns = len(df_all.columns)
    i_target = __check_column_and_target(n_columns, target_column)

    # Get the column names and the maximum value for each column
    # Here the fact that the calculation is more accurate with larger parameter is assumed.
    col_names = [x for x in df_all.columns]
    del col_names[i_target]
    col_max = []
    for col in col_names:
        col_max.append(df_all[col].max())
    t_col_name = df_all.columns[i_target]
    #print(col_names[:-1])
    #print(col_max)
    #print(df_all)

    # Group the DataFrame by groupby method
    df_all_gpb = df_all.groupby(col_names[:-1])

    # TODO:
    #   if 3D plot is required, import necessary 3D plotting modules first
    if f_plot3d:
        pass
    else:
        fig, axs = __init_fig_axs(n_columns, col_names, t_col_name)


        if n_columns == 3:
            for group in df_all_gpb.groups:
                gp_data = df_all_gpb.get_group(group)
                axs[0].plot(gp_data[col_names[-1]], gp_data[t_col_name], label="%s=%s" % (col_names[0], group[0]))

        if n_columns >= 4:
            for group in df_all_gpb.groups:
                for i in range(len(col_names)-1):
                # check the convergence of parameter col_names[i]
                # with the other parameters at the best, i.e. max
                    flag_best_other = True
                    for j in range(len(col_names)-1):
                        if j != i and group[j] != col_max[j]:
                            flag_best_other = False
                            break
                    if not flag_best_other:
                        continue

                    gp_data = df_all_gpb.get_group(group)
                    axs[i].plot(gp_data[col_names[-1]], gp_data[t_col_name], label="%s=%s" % (col_names[i], group[i]))
            # Generate the title string as the fixed parameters
            for i in range(len(col_names)-1):
                title_str_list = ['convergence w.r.t', col_names[i],'\n@ (']
                for j in range(len(col_names)-1):
                    if j != i:
                        title_str_list.append("%s = %s" % (col_names[j], col_max[j]))
                title_str_list.append(')')
                title_str = ' '.join(title_str_list)
                axs[i].set_title(title_str)

    for ax in axs:
        ax.legend(loc="upper left", shadow=True, fancybox=True)

    plt.show()

    return

# ====================================================

def Main(ArgList):

    description = ''' Visualize the data for an N-parameter convergence test. In general N=2/3, ''' 
    parser = ArgumentParser(description=description)
    parser.add_argument(dest="datafile", metavar='file', type=str, nargs=1, help="The name of file storing the data. Better in CSV/Excel form and index is not necessary.")
    #parser.add_argument("-n", dest="N", type=int, default=0, help="the number of columns in the data file (except the index). Default 0 for automatic detection.")
    parser.add_argument("-t", dest="target_column", type=int, default=0, help="the index of column (>0) which contains the quantity to converge. Default the last column.")
    parser.add_argument("--plot3d", dest="f_plot3d", action="store_true", help="Flag to use 3D plots. NOT implemented YET.")

    # initialize options as 'opts'
    opts = parser.parse_args()
    datafile = opts.datafile[0]

    pc_nd_conv_plot(datafile, opts.target_column, opts.f_plot3d)

# ==============================

if __name__ == "__main__":
    Main(sys.argv)

