# -*- coding: utf-8 -*-
'''Module that defines variables related to MYKit configuration
'''

import json
import abc
from pprint import pprint
import os
from shutil import which
import logging
import functools
 
class global_config:
    '''class for MYKit configuration
    
    MYKit uses a JSON file to manage custom variables for configuration.
    The path of the configuration JSON will first be searched by
    checking whether the environment variable ``MYKIT_CONFIG`` is defined.
    The default path for the JSON file is ``~/.mykit_config.json``.

    The built-in key-value pair has a list as the value. The first member is
    the default option value and the second is a short description.

    When setting up option in ``.mykit_config.json``, you need only to set the
    option value.
    '''
    __envVar = "MYKIT_CONFIG"
    __configPathDe = os.path.expanduser("~/.mykit_config.json")
    __options = {
        "vaspStdExec" : [which('vasp_std'), 'Path of `vasp_std` executive'],
        "mpiExec"     : [which('mpirun'), 'the MPI executable to use'],
        "numpyPrec"   : ["64", 'NumPy precision'],
        "symPrec"     : [1.0E-5, 'Symmetry tolerance for SpgLIB'],
        "vaspPAWPBE"  : [os.environ.get("VASP_POT_PBE", None), 'Path of the VASP PBE PAW directory'],
        "vaspPAWLDA"  : [os.environ.get("VASP_POT_LDA", None), 'Path of the VASP LDA PAW directory'],
        "verbWarn"    : [1, 'Verbose level of warnings'],
        "verbLog"     : [1, 'Verbose level of log information'],
        "logIndent"   : [4, 'Spaces to differentiate logs at different depth'],
    }
    __optKeys = __options.keys()

    @classmethod
    def _get_dejson_path(cls):
        '''Get the path of the JSON configuration file
        '''
        return cls.__configPathDe

    @classmethod
    def _env_var(cls):
        '''Get the name of the environment variable to set the path of custom configuration file
        '''
        return cls.__envVar

    @classmethod
    def __get_opt(cls, key, doc=False):
        assert isinstance(doc, bool)
        _v = cls.__options.get(key, None)
        if _v is not None:
            return _v[int(doc)]
        return {True: '', False: None}[doc]

    @classmethod
    def get(cls, *opts, doc=False):
        '''Get the option from the configuration.

        It supports read in multiple options and returns a tuple of the option values.
        It is recommended to get all the needed options with one ``get`` call.
        
        Args :
            *opts (str) : any number of options to extract from configuration
            doc (bool) : set True to get the description of the option, instead of value

        Returns :
            None if no opts is input, a dictionary value if opts has one argument
            and a tuple of values if opts has more than one argument.
            If ``doc`` is set to True, an empty string instead of None will be returned.
        '''
        if len(opts) == 0:
            return {True: '', False: None}[doc]
        try:
            __path = os.environ[cls._env_var()]
        except KeyError:
            __path = cls.__configPathDe
        if os.path.isfile(__path):
            logging.info("Load custrom from {}".format(__path))
            try:
                with open(__path, 'r') as _h:
                    __json = json.load(_h)
                for key in __json:
                    if key in cls.__optKeys:
                        cls.__options[key][0] = __json[key]
            except json.JSONDecodeError as _err:
                logging.warn("Fail to load custrom configuration from {}. Use default".format(__path))
                # Error in decoding the JSON config
                pass
        if len(opts) == 1:
            return cls.__get_opt(opts[0], doc=doc)
        return tuple(map(functools.partial(cls.__get_opt, doc=doc), opts))

    @classmethod
    def print_opts(cls):
        '''print out all configurable options'''
        pprint(cls.__options.keys())

