# -*- coding: utf-8 -*-
'''defines variables related to MYKit configuration
'''

import abc
import functools
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
        "vaspStdExec" : [which('vasp_std'), 'Path of `vasp_std` executive'],
        "mpiExec"     : [which('mpirun'), 'the MPI executable to use'],
        "numpyPrec"   : ["64", 'NumPy precision'],
        "symPrec"     : [1.0E-5, 'Symmetry tolerance for SpgLIB'],
        "vaspPawPbe"  : [os.environ.get("VASP_POT_PBE", None), 'Path of the VASP PBE PAW directory'],
        "vaspPawLda"  : [os.environ.get("VASP_POT_LDA", None), 'Path of the VASP LDA PAW directory'],
        "verbWarn"    : [1, 'Verbose level of warnings'],
        "verbLog"     : [1, 'Verbose level of log information'],
        "logIndent"   : [4, 'Spaces to differentiate logs at different depth'],
    }
    _optKeys = _options.keys()

    def __init__(self):
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
        # print(self._configPath)
        # print(os.path.isfile(path))
        if os.path.isfile(path):
            logging.info("Load custom from {}".format(path))
            # print("Load custom from {}".format(path))
            try:
                with open(path, 'r') as h:
                    _j = json.load(h)
                for key in _j:
                    if key in self._optKeys:
                        self._options[key][0] = _j[key]
            except json.JSONDecodeError as _err:
                # Error in decoding the JSON config
                logging.warn("Fail to load custrom configuration from {}. Use default".format(path))
                # print("Fail to load custrom configuration from {}. Use default".format(path))

    def _get_opt(self, key, doc=False):
        assert isinstance(doc, bool)
        try:
            v = self._options[key]
        except KeyError:
            raise ConfigError("Unknown option: {}".format(key))
        return v[int(doc)]

    def get(self, *opts, doc=False):
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
            return {True: '', False: None}[doc]
        if len(opts) == 1:
            return self._get_opt(opts[0], doc=doc)
        return tuple(map(functools.partial(self._get_opt, doc=doc), opts))

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
            return v[1]
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

    # @classmethod
    # def set_custom(cls, pathJson):
    #     '''Set the location of the custom configuration file
    #     '''
    #     try:
    #         assert os.path.isfile(pathJson)
    #     except AssertionError:
    #         raise ConfigError("Input config JSON not found: {}".format(pathJson))
    #     cls._configPath = pathJson
    
    # @classmethod
    # def unset_custom(cls):
    #     '''Unset the custom configuration
    #     '''
    #     cls._configPath = None
