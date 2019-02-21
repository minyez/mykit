# -*- coding: utf-8 -*-
'''Define the class for manipulating POSCAR, the VASP lattice input
'''
from mykit.core.lattice import lattice

class poscar(lattice):
    '''The class to manipulate POSCAR, the VASP lattice input file.
    '''

    def __init__(self, cell, atoms, pos, **kwargs):
        super(poscar, self).__init__(cell, atoms, pos, **kwargs)

    
    @classmethod
    def read_from_file(cls, poscarPath):
        '''Read from a existing POSCAR file.
        '''
        pass

    @classmethod
    def create_from_lattice(cls, latt):
        '''Create POSCAR from ``lattice`` instance ``latt``.
        '''
        assert isinstance(latt, lattice)
        __kw = latt.get_kwargs()
        return cls(latt.cell, latt.atoms, latt.pos, **__kw)
    