# Roadmap of `mykit`

## `core` package: Base classes


### `_control` module

Functions and classes to deal with tag and value mapping. 
They lie the basis for input conversion

  - [x] tags mapping
  - [ ] values mapping
  - [x] store all tags (keywords) mapping in a metadata file, like JSON or rST (see ASE)
  - [x] abstract base class to define the basic behavior of controller

### `unit` module

Classes that control the allowed units and their conversions inbetween

  - [x] `EnergyUnit` 
  - [x] `LengthUnit`

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
    - [x] read structure from CIF file by `PyCIFRW`

### Input controllers

Input controller classes are named after `*_control`, and are subclasses of abstract base class `prog_mapper` defined in the `_control` module.
They are metaclasses of the class of program-specific input files, 
and help map name and value of parameter tags between different simulation programs.
  - [x] `planewave_control`: plane-wave basis
  - [x] `xc_control`: (semi-)local exchange-correlation tags
  - [x] `ion_control`: ion relaxation
  - [ ] `scf_control`: electronic self-consistent loop
  - [x] `kmesh_control`: sampling in the reciprocal space
  - [ ] `lapw_control`
  - [ ] `mbpt_control`: many-body perturbation calculation control

### `kmesh` module

  - [x] kpath decoder and encoder to translate between kpath string and list of ends of path line segments

### `symmetry` module

  - wrapper functions to `spglib` to get symmetry information from `Cell` objects
    - [x] space group number and symbol
    - [x] symmetry operations
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
    - [ ] Predefined kpaths for space groups
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
    - [x] if the kvec is parsed, automatically determine if it represents a kpath, `isKpath` attribute
    - [x] eigenvalue unit argument
    - [x] if kpath is detected, compute lengths of line segments as plotting x. Otherwise use the index of kpoints.
    - [ ] export to data file
    - [x] return a `Dos` object by specifying a smearing

### `dos` module

  - `Dos`: analyze density of states (DOS)
    - [x] reading in energy grids, DOS
    - optional keyword argument
      - [x] partial (projected) DOS
      - [x] fermi level
    - [x] extract partial wave coefficients of a or one of projectors on atoms

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
  - `outcar`
    - `Outcar`
      - [x] initial structure as `Poscar` object
      - [x] return intermediate structure at any ionic step as `Poscar` object
      - [x] read atom types
      - [x] read kpoints coordinate and weights
      - [x] read number of plane waves at each kpoint
      - [x] an `Incar` object with some import tags read
      - [x] read eigenvalues and occupations at any ionic step and return a `BandStructure` object (no projection). Trim kpoints allowed.
    - [x] get OUTCAR value by key, `get_value`
  - `xml`
    - `Vaspxmlrun`
      - [x] analyze kpoints
      - [x] reading eigenvalues and occupation numbers in last calculation section
      - [x] reading partial DOS and wave function projection
      - [x] return a `BandStructure` object from the collected information
      - [x] return a `Dos` object from the collected information
      - [x] generate `Incar` object from incar section
      - [x] generate `Poscar` objects from initpos and finalpos section


## `wien2k` package: WIEN2k related

  - `constants`: set up constants for WIEN2k, such as fortran IO format, default RMT and R0
  - `utils`
    - get default RMT and R0 for elements
    - print a symmetry operation block with a rotation matrix and a translation vector
  - `inputs`
    - `In1`
      - [x] read from file
  - `struct`
    - `Struct`
      - [x] subclass of `Cell`, with extra keyword argument `rmt` and `r0`.

## `visualize` module: visualizing results
 
  - `Init`: a global wrapper for visualizers
    - [ ] call corresponding visualizer by type check
    - [x] keyword arguments to parse visualizers
  - `XmGrace`: manipulate grace files 
  - `BSVisualizer`: visualize band structure mainly with `matplotlib.pyplot`
    - [ ] agr file export
    - [x] set range for plotting
    - [x] choose spin component to plot, when two spin channels are available
    - [x] set if use VBM as energy zero (`align_vbm` keyword)
    - [ ] option for drawing DOS as well. The Dos object can be parsed directly by `dos` keyword, or set it to a float number as smearing and let `BandStructure` compute one.
    - [x] draw selected bands
    - [x] mark the kpoints symbols
    - [x] allow to draw the projections for selected components on atoms of selected bands if projections are available
      - [x] raise error if no projection is available
      - [x] drawn as dots on the band, with area representing the value of projections, see [this SO link](https://stackoverflow.com/a/14860958)
      - [x] drawn as a fill under on the band, with width representing the value of projections

## tools

  - [ ] decouple the argument parsing with the main program


## Miscs

  - [x] version string controled by [versioneer](https://github.com/warner/python-versioneer)
  - [x] Travis CI
  - [x] CircleCI
  - [x] Codacy
  - [x] Code coverage