# mykit

[![Travis CI](https://travis-ci.org/minyez/mykit.png?branch=master)](https://travis-ci.org/minyez/mykit)
[![CircleCI](https://circleci.com/gh/minyez/mykit/tree/master.svg?style=svg)](https://circleci.com/gh/minyez/mykit/tree/master)

Personal toolkit for manipulating input and output files of various electornic structure calculators for periodic systems.

*Note: Under extensive refactoring. Limited functionality and faaaaaar from the objectives.*

## Usage

Clone the repository to some local path, e.g. `/path/to/mykit/`. Install the package by

```
python setup.py install
```

and add to environment variables

```shell
export PATH="/path/to/mykit/tools/:$PATH"
```

Perform a quick test

```shell
cd /path/to/mykit/test
bash runtests.sh
```

or simply run `pytest` at `/path/to/mykit/`

Check the `doc/examples` for usages.

## Objectives

This module has its emphasis on the input files generation and tutorial purpose, by which means

- The writing of daily-use input files is faster, by implemented easy-to-use classes for the inputs, equipped with many convenient factory methods.
- Generate all inputs for a tutorial calculation with one command-line, with tags explained within as comments (if possible)
<!-- - A series of calculations can be run by statements within a python script, as in [ASE](https://wiki.fysik.dtu.dk/ase/). And it also provides an alternative to run in a bash script, with each command line a self-explained -->
- Conversion from input files of one atomic simulation program to those (or that) of another is nothing else than one-line command in terminal.

## Acknowledgement

- `Symmetry` class and `get_spg` method of `space_group` class is implemented based on python bindings of [spglib](https://atztogo.github.io/spglib/python-spglib.html) by @atztogo
- Special kpoints are retrieved from [Bilbao Crystallographic Server](http://www.cryst.ehu.es), explicitly by using the `KVEC` program. See [Aroyo2014](https://dx.doi.org/10.1107/S205327331303091X)
- version string controlled by using [versioneer](https://github.com/warner/python-versioneer)

<!-- Currently partially supported codes:

- [VASP](http://www.vasp.at/)
- [WIEN2K](http://susi.theochem.tuwien.ac.at/)
- [ABINIT](https://www.abinit.org/) -->
