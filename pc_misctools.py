#!/usr/bin/env python
# coding=utf-8

# ====================================================
#
#     File Name : pc_misctools.py
# Creation Date : 18-04-2018
#    Created By : Min-Ye Zhang
#       Contact : stevezhang@pku.edu.cn
#
# ====================================================

import sys

def Main(ArgList):
    # add from argparse import ArgumentParser at head
    description='''Miscellaneous tools for daily scientific work.'''
    
    parser = ArgumentParser(description=description)
    #group1 = parser.add_mutually_exclusive_group()
    
    parser.add_argument("-t", dest="task", type=str, help="the task to take")
    #group1.add_argument("-x", dest="exclu", type=int, default=0, help="template exclusive argument")
    parser.add_argument('--show',dest='showtasks', action="store_true", help="optional argument with undetermined number of paras")
    
    # initialize options as 'opts'
    opts = parser.parse_args()
    


# ==============================

if __name__ == "__main__":
    Main(sys.argv)

