import streusel
from streusel.gaussian_cube import *
import os

cube = os.getcwd() + '/ethane/ethane_epot.cube'

mol = Molecule(cube)
mol.get_efield()
mol.sample_efield()
print(mol.vol)


