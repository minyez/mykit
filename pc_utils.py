#!/usr/bin/env python
# coding=utf-8

# ====================================================
#     File Name : pc_utils.py
# Creation Date : 25-04-2018
#    Created By : Min-Ye Zhang
#       Contact : stevezhang@pku.edu.cn
# ====================================================

from __future__ import print_function
import os,re
import subprocess as sp

# ====================== PERFORM CALCULATION ======================
def common_run_calc_cmd(calc_cmd, fout=None, ferr=None):
    '''
    Run the calculation command by threading a subprocess calling calc_cmd
    '''

    if fout is None:
        ofile = sp.PIPE
    else:
        ofile = open(fout,'w')

    if ferr is None:
        efile = sp.PIPE
    else:
        efile = open(ferr,'w')

    p=sp.Popen(calc_cmd,stdout=ofile,stderr=efile,shell=True)
    p.wait()

    if not fout is None: ofile.close()
    if not ferr is None: efile.close()

# ====================== PRINT WARNING ======================
def common_print_warn(warn_str, func_level=0):
    '''
    Print warning with a specific level for the calling function
    '''
    print("  "*(func_level+1) + '- WARNING: ' + warn_str)


# ====================== CREATE DIRECTORY ======================

def common_io_checkdir(dirname=None,create=True):
    '''
    check if dirname exists, or create it
    return: the full path of target directory
    '''
    dirname = dirname.strip()

    if (dirname is None or dirname.strip() == ""):
        dirname = os.getcwd()
    elif (not os.path.exists(dirname)) and create:
        os.mkdir(dirname.strip())
    return dirname


def common_io_cleandir(dirname=None):
    '''
    check if dirname exists and is empty, or create it
    and make it contain no files
    return: the full path of target directory
    '''
    if (dirname is None or dirname.strip() == ""):
        dirname = os.getcwd()
    elif (not os.path.exists(dirname)):
        os.mkdir(dirname)
    elif (os.path.exists(dirname)):
        sp.call("rm -rf %s" % (dirname+"/*"),shell=True)
    return dirname

# ====================== CHECK LOCATION INFO ======================

def common_get_dirname(path='.'):
    '''
    return the directory location of path.
    Return: path,                         if path is a directory itself
            the parent directory of path, if path is a file
    '''
    abspath = os.path.abspath(path)
    if os.path.isdir(path):
        parentdir = abspath
    elif os.path.isfile(path):
        parentdir = abspath[:-len(os.path.basename(path))-1]
    else:
        parentdir = None
    return parentdir

# ====================== STRING CONVERSION ======================

def common_ss_conv(string, i, conv2, sep=None):
    '''
    Split the string and convert a single substring to a specified type.

    Parameters:
        string: str
            the string from which to convert value
        i: int
            the substring index in the list to be converted after splitting by sep
        conv2: type
            the type to which the substring will be converted
        sep: str
            the separators used to split the string.
    '''

    str_tmp = string.strip()
    if sep is not None:
        str_list = re.split(r'[%s]'%sep, str_tmp)
        print(str_list)
    else:
        str_list = str_tmp.split()

    try:
        return conv2(str_list[i])
    except ValueError:
        return conv2(float(str_list[i]))



