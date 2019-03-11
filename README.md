# mykit

Personal toolkit for manipulating input and output files of various electornic structure calculators for periodic systems.

*Note: Under extensive refactoring. Limited functionality and faaaaaar from the objectives.*

## Objectives

<!-- *This python module has no intention, and is definitely unable to compete with the repsected [ASE](https://wiki.fysik.dtu.dk/ase/) module.*-->

This module has its emphasis on the input files generation and tutorial purpose, by which means

- The writing of daily-use input files is faster, by implemented easy-to-use classes for the inputs, equipped with many convenient factory methods.
- Generate all inputs for a tutorial calculation with one command-line, with tags explained within as comments (if possible)
<!-- - A series of calculations can be run by statements within a python script, as in [ASE](https://wiki.fysik.dtu.dk/ase/). And it also provides an 
    alternative to run in a bash script, with each command line a self-explained -->
- Conversion from input files of one atomic simulation program to those (or that) of another is nothing else than one-line command in terminal.

## Acknowledgement

- `Symmetry` class is implemented based on python bindings of [spglib](https://atztogo.github.io/spglib/python-spglib.html) by @atztogo
- Special kpoints are retrived from [Bilbao Crystallographic Server](http://www.cryst.ehu.es), explicitly the `KVEC` program. Citations: [Aroyo2014](doi:10.1107/S205327331303091X) and [Tasci2012](doi:10.1051/epjconf/20122200009)


<!-- Currently partially supported codes:

- [VASP](http://www.vasp.at/)
- [WIEN2K](http://susi.theochem.tuwien.ac.at/)
- [ABINIT](https://www.abinit.org/) -->
