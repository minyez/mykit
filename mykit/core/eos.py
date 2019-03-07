#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''equation of state
'''

import numpy as np

# pylint: disable=no-member
class eos:
    '''defines various equations of state as staticmethods
    '''

    @staticmethod
    def Birch_Murnaghan(V, E0, V0, B0, Bp):
        '''Birch-Murnaghan Equation of State

        $$ E = E_0 + \frac{9}{16} $$

        Args:
            V: volume
            E0: energy at equillirium volume
            V0: equillirium volume
            B0: bulk modulus
            Bp: 
        '''
        return E0 + 9.0 / 16.0 * V0 * B0 * \
                ((np.power(V0/V, 2.0/3.0) - 1.0)**3*Bp + \
                 (np.power(V0/V, 2.0/3.0) - 1.0)**2*(6.0-4.0*np.power(V0/V, 2.0/3.0)))

