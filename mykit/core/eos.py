#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''defines various equations of state as staticmethods
'''

import numpy as np


# pylint: disable=no-member
def Birch_Murnaghan(V, E0, V0, B0, Bp):
    '''Birch-Murnaghan Equation of State

    $$ E = E_0 + \\frac{9 V_0 B_0}{16}(Bp((V_0/V)^(2/3) - 1.0)^3 + ((V_0/V)^(2/3)-1.0)^2(6-4(V_0/V)^(2/3))) $$

    Args:
        V: volume
        E0: energy at equillirium volume
        V0: equillirium volume
        B0: bulk modulus
        Bp: coefficient
    '''
    return E0 + 9.0 / 16.0 * V0 * B0 * \
        ((np.power(V0/V, 2.0/3.0) - 1.0)**3*Bp +
         (np.power(V0/V, 2.0/3.0) - 1.0)**2*(6.0-4.0*np.power(V0/V, 2.0/3.0)))
