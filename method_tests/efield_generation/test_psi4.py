from molecule import *
import os

print(os.getcwd())
parent_path = os.getcwd()
ethane_opt_xyz = '/home/he/bin/github_clones/clones4flockit/EquivariantMultipoleGNN/ethane_opt.xyz'
new_calc_path = os.getcwd() + '/test_ethane/'

mol = Molecule(parent_path, new_calc_path)
mol.init_from_xyz(ethane_opt_xyz)
mol.gen_esp_cube()

