# pylint: disable=missing-docstring
from mykit import core, vasp, wien2k

from ._version import get_versions

#from mykit import gpaw

__version__ = get_versions()['version']
del get_versions
