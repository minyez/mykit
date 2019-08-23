#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""classes and functions for manipulating plots by grace
"""
import numpy as np

class GraceDataSet:
    """object for controlling data set
    """
    def __init__(self, x, y, e=None, comment=None):
        self.x = x
        self.y = y
        assert len(np.shape(x)) == 1
        assert np.shape(x) == np.shape(y)
        self.e = e
        if e is not None:
            raise NotImplementedError
        self.comment = comment

    @property
    def type(self):
        """type of dataset, xy, xye"""
        t = 'xy'
        return t

    def __str__(self):
        ss = [f'@type {self.type}',]
        for x, y in zip(self.x, self.y):
            ss.append('%10.6f %10.6f' % (x, y))
        ss.append('&')
        return '\n'.join(ss)


class GraceGraph:
    """object for controlling graph
    """
    def __init__(self):
        pass

    @property
    def nsets(self):
        """int. Number of datasets"""
        return 0


class Grace():
    """object for the overall control of grace plots
    """

    def add_dataset_to_graph(self, igraph):
        """Add data set to graph
        """
        pass

    def add_graph(self):
        """Add a new graph to Grace
        """
        pass

    def adjust_graph(self, igraph):
        """adjust graph appearance
        """
    
    def export(self):
        """Export to agr file
        """
        pass
