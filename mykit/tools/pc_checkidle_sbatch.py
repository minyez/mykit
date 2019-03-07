#!/usr/bin/env python
# coding=utf-8

# ====================================================
#     File Name : pc_checkidle_sbatch.py
# Creation Date : 09-05-2018
#    Created By : Min-Ye Zhang
#       Contact : stevezhang@pku.edu.cn
# ====================================================

from __future__ import print_function
import subprocess as sp
import re, sys
import string


def __node_list_of_part(sinfo_line):
    '''
    Generate the list of nodes from a single sinfo line
    '''

    # list of original node name string
    node_list_ori = []
    # list of expanded node name string, which is returned
    node_list_exp = []

    comma_index = [-1]

    sinfo_line_words = sinfo_line.split()
    part = sinfo_line_words[0]
    state = sinfo_line_words[4]
    node_str = sinfo_line_words[-1]

    # split the sinfo_line
    for i in range(len(node_str)):
        if node_str[i] == ',':
            # tell if the comma in a "[]" to separate two node IDs or names
            # the node name usually starts with an alphabet
            if node_str[i+1] in string.letters:
                comma_index.append(i)
    comma_index.append(len(node_str))

    for i in range(len(comma_index)-1):
        node_list_ori.append(node_str[comma_index[i]+1:comma_index[i+1]])

    for node in node_list_ori:
        if not node.endswith(']'):
            # single node name
            node_list_exp.append(node)
        else:
            # several node IDs with the same name prefix
            node_to_expand_list = re.split(r'[\[\]]', node)[:-1]
            node_int_list = __expand_nodeid(node_to_expand_list[1])
            for node_int in node_int_list:
                node_list_exp.append("%s%02d" % (node_to_expand_list[0],node_int))

    return part, state, node_list_exp


def __expand_nodeid(nodeid_num_str):
    '''
    Expand a string of node ID numbers, such as [06,11-14,16-20] 
    and return all the ID numbers (integer) involved
    '''

    node_int_list = []
    # first separate by comma
    comma_sep = nodeid_num_str.split(',')
    for nodestr in comma_sep:
        if not '-' in nodestr:
            node_int_list.append(int(nodestr))
        else:
            startid, endid = [int(x) for x in nodestr.split('-')]
            for i in range(startid,endid+1):
                node_int_list.append(i)
    return node_int_list


def __avail_cores(scontrol_cmd, nodeid):

    scontrol_cmd_tmp = scontrol_cmd + [nodeid]
    scontrol_info = sp.check_output(scontrol_cmd_tmp)
    cpu_info = scontrol_info.split('\n')[1].split()
    
    cpu_alloc = int(cpu_info[0].split('=')[-1])
    cpu_error = int(cpu_info[1].split('=')[-1])
    cpu_total = int(cpu_info[2].split('=')[-1])
    cpu_avail = cpu_total - cpu_alloc - cpu_error

    return cpu_total, cpu_avail


def Main(ArgList):
    '''
    This simple program checks the available computational sources.
    Currently it only supports SBATCH job managing system.
    '''

    scontrol_cmd = ['scontrol','show','node']
    # need to check mix only, since idle is easy to view
    #avail_tag = ['mix','idle']
    avail_tag = ['mix']

    # the partitions you do not want to check
    ignore_part = ['GPU','KNL','TH_SHORT']
    
    sinfo = sp.check_output('sinfo')
    sinfo_list = sinfo.split('\n')[1:-1]
    
    # get all available partitions
    avail_part = []
    for info in sinfo_list:
        info_words = info.split()
        if info_words[4] in avail_tag:
            avail_part.append(info)
    
    print(" Partition State     nodeid     CPUs    Avail")
    print("---------------------------------------------")
    for sinfo_line in avail_part:
        part, state, node_in_part = __node_list_of_part(sinfo_line)
        # ignore some particular partition, e.g. GPU, KNL
        if part in ignore_part:
            continue
        for nodeid in node_in_part:
            cpu_total, cpu_avail = __avail_cores(scontrol_cmd, nodeid)
            print("%10s %5s %10s %8d %8d" % (part, state, nodeid, cpu_total, cpu_avail ))
    print("---------------------------------------------")

if __name__ == "__main__":
    Main(sys.argv)
