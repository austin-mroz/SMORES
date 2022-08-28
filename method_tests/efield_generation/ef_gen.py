import sys; sys.path.insert(0,'/home/he/bin/github_clones/STREUSEL/streusel')

from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation as R

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import psi4

import argparse
import pathlib
import tempfile
import os

# this function is adapted from code provided by Lukas Turcani in a test file
def _get_command_line_artuments():
    parser = argparse.ArgumentParser()
    parser.add_argument("xyz_file")
    return parser.parse_args()

# to avoid ase, we write our own function to read xyz files
def read_xyz(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()
    df = pd.DataFrame()
    for line in lines[2:]:
        line = list(filter(None, line.removesuffix('\n').split(' ')))
        nrow = pd.DataFrame([[str(line[0]), float(line[1]), float(line[2]), float(line[3])]])
        tempdf = pd.concat([df, nrow])
        df = pd.DataFrame(tempdf)
    df.columns = ['atom', 'x', 'y', 'z']

def translate_molecule(xyz):
    """ translates the molecule so that the center of mass is at (0,0,0) """
    def get_CenterOfMass(xyz):
        """ Returns the NG* values of the Center of Mass """
        atomic_mass = {'H':1.01, 'He':4.00, 'Li':6.94, 'Be':9.01, 'B':10.81, 'C':12.01,
                    'N':14.01, 'O':16.00, 'F':19.00, 'Ne':20.18, 'Na':22.99, 'Mg':24.31,
                    'Al':26.98, 'Si':28.09, 'P':30.97, 'S':32.07, 'Cl':35.45, 'Ar':39.95,
                    'K':39.10, 'Ca':40.08, 'Sc':44.96, 'Ti':47.87, 'V':50.94, 'Cr':52.00,
                    'Mn':54.94, 'Fe':55.85, 'Co':58.93, 'Ni':58.69, 'Cu':63.55, 'Zn':65.39,
                    'Ga':69.72, 'Ge':72.61, 'As':74.92, 'Se':78.96, 'Br':79.90, 'Kr':83.80,
                    'Rb':85.47, 'Sr':87.62, 'Y':88.91, 'Zr':91.22, 'Nb':92.91, 'Mo':95.94,
                    'Tc':98.00, 'Ru':101.07, 'Rh':102.91, 'Pd':106.42, 'Ag':107.87,
                    'Cd':112.41, 'In':114.82, 'Sn':118.71, 'Sb':121.76, 'Te':127.60,
                    'I':126.90, 'Xe':131.29, 'Cs':132.91, 'Ba':137.33, 'La':138.91,
                    'Ce':140.12, 'Pr':140.91, 'Nd':144.24, 'Pm':145.00, 'Sm':150.36,
                    'Eu':151.96, 'Gd':157.25, 'Tb':158.93, 'Dy':162.50, 'Ho':164.93,
                    'Er':167.26, 'Tm':168.93, 'Yb':173.04, 'Lu':174.97, 'Hf':178.49,
                    'Ta':180.95, 'W':183.84, 'Re':186.21, 'Os':190.23, 'Ir':192.22,
                    'Pt':195.08, 'Au':196.97, 'Hg':200.59, 'Tl':204.38, 'Pb':207.2,
                    'Bi':208.98, 'Po':209.00, 'At':210.00, 'Rn':222.00, 'Fr':223.00,
                    'Ra':226.00, 'Ac':227.00, 'Th':232.04, 'Pa':231.04, 'U':238.03,
                    'Np':237.00, 'Pu':244.00, 'Am':243.00, 'Cm':247.00, 'Bk':247.00,
                    'Cf':251.00, 'Es':252.00, 'Fm':257.00, 'Md':258.00, 'No':259.00,
                    'Lr':262.00, 'Rf':261.00, 'Db':262.00, 'Sg':266.00, 'Bh':264.00,
                    'Hs':269.00, 'Mt':268.00}
        total_mass = 0
        comx = 0
        comy = 0
        comz = 0

        for index, row in xyz.iterrows():
            atom_mass = atomic_mass[str(row['atom'])]
            comx += atom_mass * row['x']
            comy += atom_mass * row['y']
            comz += atom_mass * row['z']
            total_mass += atom_mass
        comx = comx/total_mass
        comy = comy/total_mass
        comz = comz/total_mass
        com = pd.DataFrame(data=[[comx, comy, comz]])
        com.columns = ['x', 'y', 'z']
        return com

    # obtain the center of mass for the molecule
    com = get_CenterOfMass(xyz)

    # translate the coordinates so that the center of mass is at (0,0,0)
    translated_mol = pd.DataFrame()
    for index, row in xyz.iterrows():
        nX = row['x'] - com['x'].values[0]
        nY = row['y'] - com['y'].values[0]
        nZ = row['z'] - com['z'].values[0]
        nrow = pd.DataFrame([[row['atom'], nX, nY, nZ]])
        tempdf = pd.concat([translated_mol, nrow])
        translated_mol = pd.DataFrame(tempdf)
    translated_mol.columns = ['atom', 'x', 'y', 'z']
    return translated_mol

def gen_grid4psi4(calc_path, res):
    grid_xyz_coords = []
    grid_elements = []

    for i in range(res):
        for j in range(res):
            for k in range(res):
                itrans = -5 + 0.2*i # xcnt
                jtrans = -5 + 0.2*j # ycnt 
                ktrans = -5 + 0.2*k # zcnt
                grid_xyz_coords.append([itrans, jtrans, ktrans])
                grid_elements.append([i, j, k])
    with open(f'{calc_path}grid.dat', 'w') as file:
        for xyz in grid_xyz_coords:
            for c in xyz:
                file.write(f'{c} ')
            file.write('\n')

def calc_esp_and_gen_cube_file(calc_path, trans_mol):
    psi4.set_options({'basis': 'aug-cc-pVDZ'})
    B_to_A = 0.529177249
    # define a psi4 molecule object using the atoms and coordinates from the original xyz file
    # we first extract the info from the translated molecule dataframe
    elements = list(trans_mol['atom'])
    coordinates = trans_mol[['x', 'y', 'z']].to_numpy()

    # define psi4 molecule object
    psi4_mol = psi4.core.Molecule.from_arrays(coordinates, elem=elements)
    psi4.core.set_output_file('output.dat', False)

    # extract EF & ESP grids
    E, wfn = psi4.prop('scf', molecule = psi4_mol, properties = ['GRID_ESP', 'GRID_FIELD'], return_wfn = True)

    psi4.set_options({'CUBEPROP_TASKS': ['ESP']})
    psi4.cubeprop(wfn)

def main():
    cli_args = _get_command_line_artuments()
    xyz = read_xyz(cli_args.xyz_file)
    trans_mol = translate_molecule(xyz)
    gen_grid4psi4(os.getcwd(), 51)
    
if __name__ == "__main__":
    main()