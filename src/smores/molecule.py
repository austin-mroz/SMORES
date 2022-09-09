import os
import pandas as pd 
import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem

#import sys; sys.path.insert(0,'/home/he/bin/github_clones/STREUSEL/streusel')
#from gaussian_cube import *

import psi4

import utilities as util


class Molecule:
    """
    class to handle molecule objects for featurization
    """

    def __init__(self, parent_path, calc_path):
        self.parent_path = parent_path
        self.calc_path = calc_path
        if not os.path.isdir(self.calc_path):
            os.makedirs(self.calc_path)

    def init_from_xyz(self, xyz_file_path):
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
            return df
        self.xyz = read_xyz(xyz_file_path)
        self.elements = list(self.xyz['atom'])
        self.coordinates = self.xyz[['x','y','z']].to_numpy()

    def init_from_smiles(self, smiles):

        def mk_xyz(self):
            with open(f'{self.calc_path}sol.xyz','w') as f:
                f.write(f'{len(self.elements)}\n')
                f.write(f'{self.smiles}\n')
                for element, coordinate in zip(self.elements, self.coordinates):
                    f.write(f'{element} {coordinate[0]} {coordinate[1]} {coordinate[2]}\n')
            
            return f'{self.calc_path}sol.xyz'
        self.smiles = smiles

        # Create a molecule with RDKit
        rdkit_mol = Chem.AddHs(Chem.MolFromSmiles(self.smiles))

        self.elements = [a.GetSymbol() for a in rdkit_mol.GetAtoms()]
        
        # Generate a conformation
        AllChem.EmbedMolecule(rdkit_mol)
        self.coordinates = rdkit_mol.GetConformer(0).GetPositions().astype(np.float32)
        xyz_file = mk_xyz(self)
        self.xyz = util.read_xyz(xyz_file)

    def gen_esp_cube(self, res=51, optimize=False):
        def translate_molecule(self):
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
            com = get_CenterOfMass(self.xyz)

            # translate the coordinates so that the center of mass is at (0,0,0)
            translated_mol = pd.DataFrame()
            for index, row in self.xyz.iterrows():
                nX = row['x'] - com['x'].values[0]
                nY = row['y'] - com['y'].values[0]
                nZ = row['z'] - com['z'].values[0]
                nrow = pd.DataFrame([[row['atom'], nX, nY, nZ]])
                tempdf = pd.concat([translated_mol, nrow])
                translated_mol = pd.DataFrame(tempdf)
            translated_mol.columns = ['atom', 'x', 'y', 'z']
            return translated_mol

        def gen_voxel_grid(self, res):
            # we want to make a box with one corner at (-5,5,5)A and one at (5,5,5)A with a resolution of 0.2 A
            grid_xyz_coords = []
            for i in range(res):
                for j in range(res):
                    for k in range(res):
                        itrans = -5 + 0.2*i
                        jtrans = -5 + 0.2*j 
                        ktrans = -5 + 0.2*k
                        grid_xyz_coords.append([itrans, jtrans, ktrans])
            # write the grid to a .dat file
            with open(self.calc_path + 'grid.dat', 'w') as file:
                for xyz in grid_xyz_coords:
                    for c in xyz:
                        file.write(str(c) + ' ')
                    file.write('\n')
            self.grid = self.calc_path + 'grid.dat'

        print('translating molecule to (0,0,0')
        trans_mol = translate_molecule(self)
        self.elements = list(trans_mol['atom'])
        self.coordinates = trans_mol[['x','y','z']].to_numpy()
        
        print('generating voxel grid')

        gen_voxel_grid(self, res)
        os.chdir(self.calc_path)
        print(os.getcwd())
        #import psi4
        psi4.set_options({'basis': 'aug-cc-pVDZ',
                            'CUBEPROP_TASKS': ['ESP'],
                            'CUBEPROP_FILEPATH': self.calc_path,
                            'reference': 'uhf',
                            })
        psi4.core.set_num_threads(14)

        psi4_mol = psi4.core.Molecule.from_arrays(self.coordinates, elem=self.elements)
        psi4.core.set_output_file(self.calc_path + 'output.dat', False)
        self.output = self.calc_path + 'output.dat'
        
        if optimize:
            print('optimizing!')
            psi4.optimize('PBE',molecule=psi4_mol)
        
        print('calculating ESP')
        E, wfn = psi4.prop('PBE', molecule=psi4_mol, properties=['GRID_ESP'], return_wfn=True)
        psi4.cubeprop(wfn)
        self.cubefile = self.calc_path + 'ESP.cube'
    
        os.chdir(self.parent_path)



