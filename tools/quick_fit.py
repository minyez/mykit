#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Facility to provide quick fitting of x-y data

For fitting of equation of state, use fitBMEOS instead.
"""

from io import StringIO
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from scipy.optimize import curve_fit
from scipy.stats import linregress
import matplotlib.pyplot as plt
import numpy as np

from mykit.core.utils import trim_after, XYFit

FIT_FUNCS = {
    0: ("Linear regression: y = k*x + b", 2, ("k", "b"), linregress),
    1: ("Linear of inverse x: y = a0/(x-a1) + a2", 3, ("a0", "a1", "a2"), XYFit.linear_inv_x),
    }

def _read_file_data(filename, x_col, y_col, start=0, end=-1):
    """read the data by blocks separated by empty lines

    Note that comments are also treated as empty lines

    Workflow:
    1. read the file
    2. remove all comments and empty lines before the first data block
    3. remove all comments and empty lines after the last data block
    4. find indice of empty lines inbetween data block. Comments inbetween are handles by np.loadtxt
    """
    with open(filename, 'r') as h:
        lines = [trim_after(l.strip(), "#") for l in h.readlines()]

    while lines != []:
        if lines[0] == '':
            del lines[0]
        else:
            break
    while lines != []:
        if lines[-1] == '':
            del lines[-1]
        else:
            break

    i = 0
    lid_empty = []
    consec_empty = False
    n = len(lines)
    while i < n:
        if lines[i] == '':
            if not consec_empty:
                lid_empty.append(i)
                consec_empty = True
        else:
            consec_empty = False
        i += 1
        
    n_data = 1
    xs = []
    ys = []
    if lid_empty == []:
        xs, ys = np.loadtxt(filename, unpack=True, usecols=(x_col, y_col))
        xs = xs[start:end+1]
        ys = ys[start:end+1]
    else:
        n_data = len(lid_empty) + 1
        lid_empty = [0,] + lid_empty + [n,]
        for i in range(n_data):
            s = StringIO('\n'.join(lines[lid_empty[i]:lid_empty[i+1]]))
            x, y = np.loadtxt(s, unpack=True, usecols=(x_col, y_col))
            xs.append(x[start:end])
            ys.append(y[start:end])
    return n_data, xs, ys


def _read_direct_data(d):
    """read directly typed data
    """
    assert len(d)%2 == 0
    nd = 1
    xs = np.array(d[0::2])
    ys = np.array(d[1::2])
    return nd, xs, ys


def fit_func(x, y, func_id, p0=None, maxfev=5000):
    """Fit
    """
    n = len(x)
    assert n == len(y)
    R2 = 0.0
    f = FIT_FUNCS[func_id][-1]
    p = []
    # linear regression use linregress
    if func_id == 0:
        p = f(x, y)
        p = p[:2]
        R2 = 1 - np.sum((p[0]*x+p[1]-y)**2)/np.sum((y-np.mean(y))**2)
    else:
        p, _ = curve_fit(f, x, y, p0=p0, maxfev=maxfev)
        R2 = 1 - np.sum((f(x, *p)-y)**2)/np.sum((y-np.mean(y))**2)
    return p, R2


def plot_fitted_func(ax, x, y, func_id, p, **kwargs):
    """plot the
    """
    _x = np.linspace(min(x), max(x), 100)
    f = FIT_FUNCS[func_id][-1]
    if func_id == 0:
        _y = p[0] * _x + p[1]
    else:
        _y = f(_x, *p)
    l = ax.plot(x, y, ls='none', marker='o', **kwargs)
    ax.plot(_x, _y, color=l[0].get_color())


def quick_fit():
    """the main stream
    """
    funcs_doc = '\n'.join(["%s. %s (%s params)" % (k, v[0], v[1]) for k, v in FIT_FUNCS.items()])
    parser = ArgumentParser(description=__doc__ + "\n" + funcs_doc, \
            formatter_class=RawDescriptionHelpFormatter)
    data_src = parser.add_mutually_exclusive_group(required=True)
    data_src.add_argument("-f", default=None, type=str, help="name of file to extract data")
    data_src.add_argument("-d", default=None, type=float, nargs="+", \
            help="direct mode for data input, (x1,y1,x2,y2...)")
    parser.add_argument("-x", dest="ix", type=int, default=0, help="<ix=0> column index of x data")
    parser.add_argument("-y", dest="iy", type=int, default=1, help="<iy=1> column index of y data")
    parser.add_argument("-s", dest="start", type=int, default=0, help="<s=0> index of start y data")
    parser.add_argument("-e", dest="end", type=int, default=int(1e9), \
            help="<e=1e9> index of end y data (inclusive)")
    parser.add_argument("-t", dest="func_id", choices=FIT_FUNCS.keys(), type=int, default=0, \
            help="index of fitting function")
    parser.add_argument("-p", dest="p0", nargs="+", help="initial guess for the parameters")
    parser.add_argument("--pp", dest="printp", nargs="+", type=int, \
            help="indices of paramters to print")
    parser.add_argument("--plot", action="store_true", help="plot the fitted curve")
    parser.add_argument("-D", dest="debug", action="store_true", help="debug mode")

    args = parser.parse_args()
    if args.f is not None:
        nd, xs, ys = _read_file_data(args.f, args.ix, args.iy, args.start, args.end)
    elif args.d is not None:
        nd, xs, ys = _read_direct_data(args.d)
    else:
        raise IOError

    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    if args.debug:
        print("# data sets: ", nd)
    p0 = args.p0
    if p0 is not None:
        p0 = [float(x) for x in p0]

    if nd == 1:
        if args.debug:
            print(xs)
            print(ys)
        p, R2 = fit_func(xs, ys, args.func_id, p0=p0)
        if args.plot:
            plot_fitted_func(ax, xs, ys, args.func_id, p)
        if args.printp is None:
            print("#points:", len(xs), ". params:", *p, " (R2=%.5f)" % R2)
        else:
            print("#points:", len(xs), ". params:", *[p[x] for x in args.printp])
    else:
        for i, (x, y) in enumerate(zip(xs, ys)):
            if args.debug:
                print(x)
                print(y)
            p, R2 = fit_func(x, y, args.func_id, p0=p0)
            if args.plot:
                plot_fitted_func(ax, x, y, args.func_id, p, label="DS %s" % i)
            if args.printp is None:
                print("#points:", len(xs), ". params:", *p, " (R2=%.5f)" % R2)
            else:
                print("#points:", len(xs), ". params:", *[p[x] for x in args.printp])

    fig.tight_layout()
    if args.plot:
        ax.legend(fontsize=16)
        plt.show()

if __name__ == "__main__":
    quick_fit()

