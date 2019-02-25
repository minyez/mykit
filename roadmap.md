# Roadmap of `MYKit`

## `core`: Base classes

- `lattice` base class: creation of periodic systems
  - [ ] Factory methods for common crystals
- Controllers: map name and value of parameter tags between different simulation programs
  - [ ] `control_map` as a base for `*_control` classes dealing with tag and value mapping
    - [x] tags mapping
    - [ ] values mapping
    - [ ] store all tags (keywords) mapping in a metadata file, like JSON or rST (see ASE)
  - [x] `planewave_control` for tags determining plane-wave basis
  - [x] `xc_control` for local exchange-correlation tags
  - [ ] `ion_control` for ion-relaxing tags
  - [ ] `scf_control` for tags controlling electronic self-consistent loop
  - [ ] `kpoint_control`
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
