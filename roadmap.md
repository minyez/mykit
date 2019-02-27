# Roadmap of `mykit`

## `core`: Base classes


### `_control` module

Functions for `*_control` classes to deal with tag and value mapping
  - [x] tags mapping
  - [ ] values mapping
  - [ ] store all tags (keywords) mapping in a metadata file, like JSON or rST (see ASE)

### `lattice` module

- `lattice` base class: creation of periodic systems
  - [ ] Factory methods for common crystals

### Input controllers

Map name and value of parameter tags between different simulation programs
  - [x] `planewave_control`: plane-wave basis
  - [x] `xc_control`: (semi-)local exchange-correlation tags
  - [ ] `ion_control`: ion relaxation
  - [ ] `scf_control`: electronic self-consistent loop
  - [ ] `kpoint_control`: sampling in the reciprocal space
  - [ ] `lapw_control`
  - [ ] `mbpt_control` for many-body perturbation calculation tags


## `vasp`: VASP related

- `incar`
  - [x] Read, print and write 
  - [ ] Tag explanation
  - [ ] Move tags into a metadata file
- `poscar` 
  - [x] Read, print and write
- `potcar`
- `potcar_search` for easily searching POTCAR
- `wavecar`
- `outcar`
- `chgcar`
- `procar`
- `contcar`
