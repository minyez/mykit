# Roadmap of `mykit`

## `core`: Base classes


### `_control` module

Functions and classes to deal with tag and value mapping. 
They lie the basis for input conversion
  - [x] tags mapping
  - [ ] values mapping
  - [x] store all tags (keywords) mapping in a metadata file, like JSON or rST (see ASE)

### `lattice` module

- `lattice` base class: creation of periodic systems
  - [ ] Factory methods for common crystals

### Input controllers

Input controller classes are named as `*_control`, and are subclasses of abstract base
classes defined in the `_control` module.
They are metaclasses of the class of program-specific input files, 
and help map name and value of parameter tags between different simulation programs.
  - [x] `planewave_control`: plane-wave basis
  - [x] `xc_control`: (semi-)local exchange-correlation tags
  - [x] `ion_control`: ion relaxation
  - [ ] `scf_control`: electronic self-consistent loop
  - [ ] `kmesh_control`: sampling in the reciprocal space
  - [ ] `lapw_control`
  - [ ] `mbpt_control`: many-body perturbation calculation control


## `vasp`: VASP related

- `incar`
  - [x] Read, print and write 
  - [ ] Tag explanation
  - [x] Move tags into a metadata file
- `poscar` 
  - [x] Read, print and write
- `potcar`
- `potcar_search` for easily searching POTCAR
- `wavecar`
- `outcar`
- `chgcar`
- `procar`
- `contcar`
