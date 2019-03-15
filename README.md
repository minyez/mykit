# mykit

[![Travis CI](https://travis-ci.org/minyez/mykit.svg?branch=master)](https://travis-ci.org/minyez/mykit)
[![CircleCI](https://circleci.com/gh/minyez/mykit/tree/master.svg?style=svg)](https://circleci.com/gh/minyez/mykit/tree/master)
[![MIT License](https://img.shields.io/badge/License-MIT-blue.svg)][https://github.com/minyez/mykit/blob/master/LICENSE]

Personal toolkit for manipulating input and output files of various electornic structure calculators for periodic systems.

*Note: Under extensive refactoring. Limited functionality and faaaaaar from the objectives.*

## Requirements

- Python >= 3.6.6
- NumPy >= 1.15.0
- Spglib >= 1.12

## Objectives

This module has its emphasis on the input files generation and tutorial purpose, by which means

- The writing of daily-use input files is faster, by implemented easy-to-use classes for the inputs, equipped with many convenient factory methods.
- Generate all inputs for a tutorial calculation with one command-line, with tags explained within as comments (if possible)
<!-- - A series of calculations can be run by statements within a python script, as in [ASE](https://wiki.fysik.dtu.dk/ase/). And it also provides an alternative to run in a bash script, with each command line a self-explained -->
- Conversion from input files of one atomic simulation program to those (or that) of another is nothing else than one-line command in terminal.

## Install

Clone the repository to some local path, e.g. `/path/to/mykit/`. 
`cd` too `/path/to/mykit/`, install the package by

```
python setup.py install
```

and add to environment variables

```shell
export PATH="/path/to/mykit/tools/:$PATH"
```

Perform a quick test to check the installation

```shell
cd /path/to/mykit/test
bash run_unittest.sh
```

or simply run `pytest` at `/path/to/mykit/`

## Usage

After installation, the jupyter notebooks in `doc/examples` may be checked for usage.

## Acknowledgement

- `Symmetry` class is implemented based on python bindings of [spglib](https://atztogo.github.io/spglib/python-spglib.html) by @atztogo
- Special kpoints are retrieved from [Bilbao Crystallographic Server](http://www.cryst.ehu.es), explicitly by using the `KVEC` program. See [Aroyo2014](https://dx.doi.org/10.1107/S205327331303091X)
- Version string controlled by using [versioneer](https://github.com/warner/python-versioneer)

<!-- Currently partially supported codes:

- [VASP](http://www.vasp.at/)
- [WIEN2K](http://susi.theochem.tuwien.ac.at/)
- [ABINIT](https://www.abinit.org/) -->
