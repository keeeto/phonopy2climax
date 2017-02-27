# phonopy2climax
## Purpose
This script takes a `phonopy` file, `mesh.yaml` and extracts information on the eigenvectors, writing them out in a format that is compatible with the a-climax code, for calculating inelastic neutron scattering spectra.

## Requirements
The code relies on a number of external python modules
* yaml
* ase
* numpy

## Use
To use the script, you must first have a `mesh.yaml` file with the eigenvectors printed out - this is achieved by setting the `EIGENVECTORS=.TRUE.` tag in your `mesh.conf` file. You must also have the original coordinates used for building your phonon calculation in `VASP` `POSCAR` format.

The script is then run as `python phonopy2climax.py`
