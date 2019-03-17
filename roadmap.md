# Roadmap of `mykit`

## `core`: Base classes


### `_control` module

Functions and classes to deal with tag and value mapping. 
They lie the basis for input conversion
  - [x] tags mapping
  - [ ] values mapping
  - [x] store all tags (keywords) mapping in a metadata file, like JSON or rST (see ASE)
  - [x] abstract base class to define the basic behavior of controller

### `cell` module

- `Cell` base class: creation of periodic systems
  - Factory methods for common crystals
    - [x] cubic cells (cP, cI, cF)
    - [x] orthorhombic cells (oP, oI, oF)
    - [ ] hexagonal and trigonal cells
    - [x] zincblende
    - [x] wurtzite
    - [x] rutile
    - [x] anatase
    - [x] diamond
    - [x] perovskite
  - [x] allow adding atoms 'on the fly'

### Input controllers

Input controller classes are named as `*_control`, and are subclasses of abstract base class `prog_mapper` defined in the `_control` module.
They are metaclasses of the class of program-specific input files, 
and help map name and value of parameter tags between different simulation programs.
- [x] `planewave_control`: plane-wave basis
- [x] `xc_control`: (semi-)local exchange-correlation tags
- [x] `ion_control`: ion relaxation
- [ ] `scf_control`: electronic self-consistent loop
- [x] `kmesh_control`: sampling in the reciprocal space
  - [x] decide the names of necessary tags
- [ ] `lapw_control`
- [ ] `mbpt_control`: many-body perturbation calculation control

### `kmesh` module

- [x] kpath decoder and encoder to translate between kpath string and list of ends of path line segments

### `symmetry` module

- `Symmetry` class for symmetry information of lattice cell, with help of [spglib](https://atztogo.github.io/spglib/python-spglib.html)
  - [x] irreducible kpoints
  - [x] return primitive cell
  - [x] return standardized cell
- `special_kpoints`
  - [x] implement special kpoints in reciprocal lattice vector of primitive cell for each space group from Bilbao server
  - [ ] k-point path from [Setyawan and Curtarolo](https://doi.org/10.1016/j.commatsci.2010.05.010)
    - [ ] create special kpoints from it, alternative to Bilbao
  - [x] converting primitive to conventional coordinate
    - [x] transformation matrix from primitive to conventional
  - [x] converting kpath string to symbols and coordinates of corresponding points on path line segments
  - [ ] Predefined kpaths for space groups :wrench:
  - [x] Factory methods for kpath of special kpoints

## `vasp` package: VASP related

- `incar`
  - [x] Read, print and write 
  - [ ] Tag explanation
  - [x] Move tags into a metadata file
- `poscar` 
  - [x] Read, print and write
- `potcar`
  - [x] `potcar_search` for easily searching POTCAR
- `kpoints`
  - [x] Decide the way to parse tags: use `kmesh_control` as attribute, value check at initialization
  - [x] `__repr__` and `__str__` magics
- `wavecar`
- `chgcar`
- `procar`
- `contcar`


## tools

- [ ] decouple the argument parsing with the main program


## Miscs

- [x] version string controled by [versioneer](https://github.com/warner/python-versioneer)
- [x] Travis CI
- [x] CircleCI
- [x] Codacy
- [x] Code coverage