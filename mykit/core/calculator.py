# -*- coding: utf-8 -*-
'''
define calculator metaclass
'''

import subprocess as sp


class calculator:
    '''Calculator metaclass

    TODO:
        add output and error file argument
    '''

    runcmd = []

    def run(self, background=False):
        '''
        Run the calculation with defined command
        '''
        calc = sp.Popen(self.runcmd)
        if not background:
            calc.wait()


