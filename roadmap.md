# Roadmap of `mykit`

## `core` package: Base classes


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
    - [x] zincblende and wurtzite
    - [x] rutile and anatase
    - [x] pyrite and marcasite
    - [x] diamond
    - [x] perovskite
  - [x] allow adding atoms 'on the fly'

### Input controllers

Input controller classes are named after `*_control`, and are subclasses of abstract base class `prog_mapper` defined in the `_control` module.
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
- `SpaceGroup`
  - [x] space group names (to cross-check with Spglib)
  - [x] tranformation matrix of kpoints coordinates from primitive to conventional cell
  - [ ] grouping of space groups according to point groups
- `SpecialKpoints`
  - [x] implement special kpoints in reciprocal lattice vector of primitive cell for each space group from Bilbao server
  - [ ] k-point paths from [Setyawan and Curtarolo](https://doi.org/10.1016/j.commatsci.2010.05.010)
    - [ ] create special kpoints from it, alternative to Bilbao
  - [x] converting primitive to conventional coordinate
    - [x] transformation matrix from primitive to conventional
  - [x] converting kpath string to symbols and coordinates of corresponding points on path line segments
  - [ ] Predefined kpaths for space groups :wrench:
  - [x] Factory methods for kpaths of special kpoints

### `bandstructure` module

- `BandStructure`: analyze band structure
  - [x] reading in eigenvalues, occupation numbers and kpoints weight
  - optional keyword arguments
    - [x] partial waves (dict)
    - [x] Fermi level `efermi` (float)
    - [x] kpoint vectors `kvec` for evaluating length of kpath
  - [x] find values and indices of VBM/CBM
    - [x] at each spin-kpoint channel
    - [x] for each spin
    - [x] overall
  - [x] compute direct band gap, fundamental gap and k-averaged gap
  - [x] extract partial wave coefficients of a or one of projectors on atoms
  - [x] compute effective gap with projections
  - [ ] if the kvec is parsed, automatically determine if it represents a kpath
  - [ ] return a `Dos` object

### `dos` module

- `Dos`: analyze density of states (DOS)
  - [ ] reading in energy grids, DOS
  - optional keyword argument
    - [ ] partial (projected) DOS
    - [ ] fermi level

## `vasp` package: VASP related

- `incar`
  - `Incar`
    - [x] Read, print and write 
    - [ ] Tag explanation
    - [x] Move tags into a metadata file
- `poscar` 
  - `Poscar`
    - [x] Read, print and write
- `potcar`
  - [x] `PotcarSearch` for easily searching POTCAR
- `kpoints`
  - `Kpoints`
    - [x] Decide the way to parse tags: use `kmesh_control` as attribute, value check at initialization
    - [x] `__repr__` and `__str__` magics
- `wavecar`
- `xml`
  - `Vaspxmlrun`
    - [x] analyze kpoints
    - [x] reading eigenvalues and occupation numbers in last calculation section
    - [x] reading partial DOS and wave function projection
    - [ ] return a `BandStructure` object from the collected information
    - [x] generate `Incar` object from incar section
    - [x] generate `Poscar` objects from initpos and finalpos section


## `wien2k` package: WIEN2k related

- `inputs`
  - `In1`
    - [x] read from file
- `struct`
  - `Struct`

## `visualize` package: visualizing results

- Matplotlib
- XmGrace

## tools

- [ ] decouple the argument parsing with the main program


## Miscs

- [x] version string controled by [versioneer](https://github.com/warner/python-versioneer)
- [x] Travis CI
- [x] CircleCI
- [x] Codacy
- [x] Code coverage