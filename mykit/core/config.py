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
    _envVar = "MYKIT_CONFIG"
    _configPathDe = os.path.expanduser("~/.mykit_config.json")
    _options = {
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
    _optKeys = _options.keys()

    @classmethod
    def _get_dejson_path(cls):
        '''Get the path of the JSON configuration file
        '''
        return cls._configPathDe

    @classmethod
    def env_var(cls):
        '''Get the name of the environment variable to set the path of custom configuration file
        '''
        return cls._envVar

    @classmethod
    def _get_opt(cls, key, doc=False):
        assert isinstance(doc, bool)
        v = cls._options.get(key, None)
        if v is not None:
            return v[int(doc)]
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
            path = os.environ[cls.env_var()]
        except KeyError:
            path = cls._configPathDe
        if os.path.isfile(path):
            logging.info("Load custrom from {}".format(path))
            try:
                with open(path, 'r') as h:
                    _j = json.load(h)
                for key in _j:
                    if key in cls._optKeys:
                        cls._options[key][0] = _j[key]
            except json.JSONDecodeError as _err:
                logging.warn("Fail to load custrom configuration from {}. Use default".format(path))
                # Error in decoding the JSON config
                pass
        if len(opts) == 1:
            return cls._get_opt(opts[0], doc=doc)
        return tuple(map(functools.partial(cls._get_opt, doc=doc), opts))

    @classmethod
    def print_opts(cls):
        '''print out all configurable options'''
        pprint(cls._options.keys())

