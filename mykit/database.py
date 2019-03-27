# coding=utf-8
'''defines classes for building models from your own database
'''
import os

from mykit.core.config import config
from mykit.core.log import Verbose


class DatabaseError(Exception):
    pass


class _db:
    '''base class for database
    '''
    pass


class Structure(_db, Verbose):

    db = config.get("dbStructure")
    if db is None:
        raise ValueError("No structure database is found")
    if not os.path.isdir(db):
        raise ValueError("No structure database is found")

    raise NotImplementedError



del config
