from horton import *
import numpy as np
from mypylib import *
np.set_printoptions(linewidth=150)


## Load the coordinates from file.
mol = IOData.from_file('../../docs/h2o.xyz')

## Create a Gaussian basis set
obs = get_gobasis(mol.coordinates, mol.numbers, 'STO-3G')
nbf = obs.nbasis

## Compute Gaussian integrals
olp = obs.compute_overlap()
kin = obs.compute_kinetic()
na = obs.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers)
er = obs.compute_electron_repulsion() # These integrals are given in CHEMIST'S notation: (11|22)

# Convert the computed integrals to PHYSICIST'S notation <12|12>
eri = np.swapaxes(er, 1, 2)


## Save the matrices to a .data-file
to_file(olp, 'overlap.data')
to_file(kin, 'kinetic.data')
to_file(na, 'nuclear.data')
to_file(eri, 'two_electron.data')
