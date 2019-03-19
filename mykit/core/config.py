# -*- coding: utf-8 -*-
'''defines variables related to MYKit configuration
'''

import abc
import json
import logging
import os
from pprint import pprint
from shutil import which


class ConfigError(Exception):
    pass


class global_config:
    '''class for mykit configuration
    
    mykit uses a JSON file to manage custom variables for configuration.
    The default path for the JSON file is ``~/.mykit_config.json``.
    Custom configuration file can be specified by "MYKIT_CONFIG" environment variable

    The built-in key-value pair has a list as the value. The first member is
    the default option value and the second is a short description.

    When setting up option in ``.mykit_config.json``, you need only to set the
    option value.
    '''
    _envVar = "MYKIT_CONFIG"
    _configPathDe = os.path.expanduser("~/.mykit_config.json")
    _options = {
        "vaspStdExec" : 'Path of `vasp_std` executive',
        "mpiExec"     : 'the MPI executable to use',
        "numpyPrec"   : 'NumPy precision',
        "symPrec"     : 'Symmetry tolerance for SpgLIB',
        "vaspPawPbe"  : 'Path of the VASP PBE PAW directory',
        "vaspPawLda"  : 'Path of the VASP LDA PAW directory',
        "verbWarn"    : 'Verbose level of warnings',
        "verbLog"     : 'Verbose level of log information',
        "logIndent"   : 'Spaces to differentiate logs at different depth',
    }
    _optKeys = _options.keys()

    def __init__(self):
        self._options = {
            "vaspStdExec" : which('vasp_std'),
            "mpiExec"     : which('mpirun'),
            "numpyPrec"   : "64",
            "symPrec"     : 1.0E-5, 
            "vaspPawPbe"  : os.environ.get("VASP_POT_PBE", None),
            "vaspPawLda"  : os.environ.get("VASP_POT_LDA", None),
            "verbWarn"    : 1,
            "verbLog"     : 1,
            "logIndent"   : 4,
        }
        # try manual input config file, otherwise search for the custom file
        try:
            path = os.environ[self._envVar]
            if not os.path.isfile(path):
                raise FileNotFoundError
        except KeyError:
            path = self._configPathDe
        except FileNotFoundError:
            raise ConfigError("MYKIT_CONFIG file not exist: {}".format(path))
        self._configPath = path
        if os.path.isfile(path):
            logging.info("Load custom from {}".format(path))
            # print("Load custom from {}".format(path))
            try:
                with open(path, 'r') as h:
                    _j = json.load(h)
                for key in _j:
                    if key in self._optKeys:
                        self._options[key] = _j[key]
            except json.JSONDecodeError as _err:
                # Error in decoding the JSON config
                logging.warn("Fail to load custrom configuration from {}. Use default".format(path))
                # print("Fail to load custrom configuration from {}. Use default".format(path))

    def _get_opt(self, key):
        try:
            v = self._options[key]
        except KeyError:
            raise ConfigError("Unknown option: {}".format(key))
        return v

    def get(self, *opts):
        '''Get the option value from the configuration

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
            return None
        if len(opts) == 1:
            return self._get_opt(opts[0])
        return tuple(map(self._get_opt, opts))

    @property
    def configFrom(self):
        return self._configPath

    @classmethod
    def _get_dejson_path(cls):
        '''Get the path of the JSON configuration file
        '''
        return cls._configPathDe

    @classmethod
    def env_var(cls):
        return cls._envVar

    @classmethod
    def _get_doc_one(cls, key):
        v = cls._options.get(key, None)
        if v is not None:
            return v
        return ''

    @classmethod
    def get_doc(cls, *opts):
        '''Get the doc string of options

        Args:
            opts (str): the option to extract docstring. 
                If empty, all available options will be printed.

        '''
        _ret = None
        if len(opts) == 0:
            _ret = tuple(map(cls._get_doc_one, cls._optKeys))
        else:
            _ret = tuple(map(cls._get_doc_one, opts))
        if len(opts) == 1:
            _ret = _ret[0]
        return _ret

    @classmethod
    def print_opts(cls):
        '''print out all configurable options'''
        pprint(cls._options.items())


# used for static configuration
config = global_config()
