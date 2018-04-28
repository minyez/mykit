#!/usr/bin/env python
# coding=utf-8

# ====================================================
#     File Name : pc_misctools.py
# Creation Date : 18-04-2018
#    Created By : Min-Ye Zhang
#       Contact : stevezhang@pku.edu.cn
# ====================================================

from __future__ import print_function, absolute_import
import sys
import subprocess as sp
from argparse import ArgumentParser

def __gen_unit_alias_dictionary(unit_alias_dict, std_unit_name):

    '''
    Generate the standard-unit to alias dictionary.
    An empty dictionary is generated if std_unit_name is turned on
    '''
    alias_unit_dict = {}
    if not std_unit_name:
        for key in unit_alias_dict.keys():
            for unit in unit_alias_dict[key]:
                alias_unit_dict.update({unit: key})
    return alias_unit_dict


def __gen_std_conv_str(conv_str, alias_unit_dict, std_unit_name):

    if len(conv_str) == 1:
        flag_from_to_units = False
        units_input = conv_str[0].split('2')
        from_unit = units_input[0].lower()
        to_unit   = units_input[1].lower()
    elif len(conv_str) == 2:
        from_unit = conv_str[0].lower()
        to_unit = conv_str[1].lower()
    else:
        raise IOError("A2B does not accept arguments larger than 2.")

    if not std_unit_name:
        std_from_unit = alias_unit_dict[from_unit]
        std_to_unit   = alias_unit_dict[to_unit]
        std_conv_str  = std_from_unit + '2' + std_to_unit
    else:
        std_conv_str  = conv_str[0]
        std_from_unit = from_unit 
        std_to_unit   = to_unit
    return flag_from_to_units, std_conv_str, std_from_unit, std_to_unit

# ====================================================

def comm_conv_length(conv_str, value=1, use_gnu=False, std_unit_name=False, print_result=False, debug=False):

    '''
    Length unit conversion. Use std_unit_name to avoid dictionary converting
    '''
    # define the possible alias of units
    atomiclength_alias_list = ['atomiclength','bohr','au', 'a0']
    angstrom_alias_list     = ['angstrom', 'ang', 'a']
    m_alias_list            = ['m','meter']

    alias_lists = [atomiclength_alias_list, angstrom_alias_list, m_alias_list]

    length_unit_alias_dict = {}
    for alias_list in alias_lists:
        length_unit_alias_dict.update({alias_list[0]: alias_list})

    # define the conversion coefficient
    # if GNU-units is not used
    if not use_gnu:
        from pc_constants import abohr
        conv_coeffs_list = {
                        'angstrom2atomiclength' : 1.0E-10/abohr , \
                        'atomiclength2angstrom' : abohr*1.0E10  , \
                        'm2angstrom'            : 1.0E10        , \
                        'angstrom2m'            : 1.0E-10       , \
                        'atomiclength2m'        : abohr         , \
                        'm2atomiclength'        : 1.0/abohr     , \
                           }
    if debug:
        print(length_alias_unit_dict)

    # Generate alias-to-std-unit dictionary
    length_alias_unit_dict = __gen_unit_alias_dictionary(length_unit_alias_dict, std_unit_name)
    flag_from_to_units, std_conv_str, std_from_unit, std_to_unit \
              = __gen_std_conv_str(conv_str, length_alias_unit_dict, std_unit_name)
  
    return __return_converted_value(value, std_conv_str, conv_coeffs_list, \
                                    std_from_unit, std_to_unit, use_gnu, print_result, debug)

# ====================================================

def comm_conv_energy(conv_str, value=1, use_gnu=False, std_unit_name=False, print_result=False, debug=False):

    '''
    Energy unit conversion. Use std_unit_name to avoid dictionary converting
    '''
    # define the possible alias of units
    ev_alias_list   = ['ev']
    ha_alias_list   = ['hatree', 'ha', 'au']
    ry_alias_list   = ['rydberg', 'ry', 'ryd']
    kcal_alias_list = ['kcal']
    j_alias_list    = ['j']
    kj_alias_list   = ['kj']

    alias_lists = [ev_alias_list, ha_alias_list, ry_alias_list, kcal_alias_list, \
                   j_alias_list, kj_alias_list]

    energy_unit_alias_dict = {}
    for alias_list in alias_lists:
        energy_unit_alias_dict.update({alias_list[0]: alias_list})

    # define the conversion coefficient
    # if GNU-units is not used
    if not use_gnu:
        from pc_constants import eV2Ha, Ha2eV, eV2Ry, Ry2eV, eV2J, eV2kJ
        conv_coeffs_list = {
                        'ev2hatree'      : eV2Ha     , \
                        'hatree2ev'      : Ha2eV     , \
                        'ev2rydberg'     : eV2Ry     , \
                        'rydberg2ev'     : Ry2eV     , \
                        'hatree2rydberg' : 2.0       , \
                        'rydberg2hatree' : 0.5       , \
                        'ev2j'           : eV2J      , \
                        'j2ev'           : 1.0/eV2J  , \
                        'ev2kj'          : eV2kJ     , \
                        'kj2ev'          : 1.0/eV2kJ , \
                           }

    # Generate alias-to-std-unit dictionary
    energy_alias_unit_dict = __gen_unit_alias_dictionary(energy_unit_alias_dict, std_unit_name)
    flag_from_to_units, std_conv_str, std_from_unit, std_to_unit \
              = __gen_std_conv_str(conv_str, energy_alias_unit_dict, std_unit_name)
  
    return __return_converted_value(value, std_conv_str, conv_coeffs_list, \
                                    std_from_unit, std_to_unit, use_gnu, print_result, debug)

# ====================================================

def __return_converted_value(value, std_conv_str, conv_coeffs_list, \
                             std_from_unit, std_to_unit, use_gnu=False, print_result=False, debug=False):

    value_converted = 1.0
    if not use_gnu:
        value_converted = value * conv_coeffs_list[std_conv_str]
    else:
        out_gnu_units = sp.check_output(['units','%s %s' % (value, std_from_unit), '%s' % std_to_unit])
        if debug:
            print(out_gnu_units.split())
            print(flag_from_to_units, std_conv_str, std_from_unit, std_to_unit)
        value_converted = float(out_gnu_units.split()[1])

    if print_result:
        print(" %s %s = %e %s" % \
              (value, std_from_unit, value_converted, std_to_unit) )
    
    return value_converted

# ====================================================

def Main(ArgList):
    description='''Convert units by using scipy.constants as default'''
    
    parser = ArgumentParser(description=description)
    
    parser.add_argument('conv_str', metavar='A2B', nargs='+', help="A2B, a single string for conversion, e.g. 'au2a', or two units name, e.g ['au','a'], to convert Bohr to Angstrom. Capitalization-insensitive.")
    parser.add_argument("-t", dest='task', type=str, default='', help="the type of unit conversion to use. Use GNU units if absent.")
    parser.add_argument("-v", dest='value', type=float, default=1.0, help="Value to convert. Defaul 1")
    #parser.add_argument("--ft", dest='unit', nargs=2, type=str, help="from-unit and to-unit, e.g. au a")
    parser.add_argument("--gnu", action='store_true', help="flag to use GNU-units conversion for complex units")
    parser.add_argument("--std", dest='std_unit_name', action='store_true', help="flag to use GNU-units conversion for complex units")
    parser.add_argument("-D", dest='debug', action='store_true', help="Flag for Debug mode")

    # initialize options as 'opts'
    opts = parser.parse_args()

    '''
    Set if the GNU-units will be used
    '''
    if opts.gnu or opts.task=='':
        use_gnu = True
    else:
        use_gnu = False

    conv_str = opts.conv_str

    if opts.task != '':
        '''
        use the initial letter to determine the type of unit conversion
        '''
        task = opts.task.lower()[0]

    if task == 'l':
        comm_conv_length(conv_str, opts.value, use_gnu, print_result=True, std_unit_name=opts.std_unit_name, debug=opts.debug)
    elif task == 'e':
        comm_conv_energy(conv_str, opts.value, use_gnu, print_result=True, std_unit_name=opts.std_unit_name, debug=opts.debug)

    return
# ==============================

if __name__ == "__main__":
    Main(sys.argv)

