import dbstep.Dbstep as db
import os

log_file = os.getcwd() + '/ethane/ethane.log'

mol = db.dbstep(log_file, atom1=2, atom2=5, commandline=True, verbose=True, sterimol=True, measure='classic')

print(mol.L)
print(mol.Bmin)
print(mol.Bmax)

cube_file = os.getcwd() + '/ethane/ethane_epot.cube'
mol = db.dbstep(cube_file, atom1=2, atom2=5, commandline=True, verbose=True, sterimol=True, measure='grid', surface='vdw')

print(mol.L)
print(mol.Bmin)
print(mol.Bmax)


