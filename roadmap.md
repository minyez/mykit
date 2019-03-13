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
  - Factory methods for common crystals, check [this site](http://www.kayelaby.npl.co.uk/chemistry/3_7/3_7_7.html)
    - [x] BCC, FCC
    - [x] zincblende
    - [x] wurtzite
    - [x] rutile
    - [x] diamond
    - [x] perovskite

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

- [x] kpath decoder and encoder to translate between kpath string and list of path segments

### `symmetry` module

- `Symmetry` class for symmetry information of lattice cell, with help of [spglib](https://atztogo.github.io/spglib/python-spglib.html)
  - [x] irreducible kpoints
  - [x] return primitive cell
  - [x] return standardized cell
- `special_kpoints`
  - [x] implement special kpoints in reciprocal lattice vector of primitive cell for each space group from Bilbao server
  - [ ] converting primitive to conventional coordinate
    - [x] transformation matrix from primitive to conventional
  - [ ] Factory methods for kpath of special kpoints

## `vasp`: VASP related

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
- `wavecar`
- `chgcar`
- `procar`
- `contcar`
