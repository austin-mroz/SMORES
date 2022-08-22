import streusel
from streusel.gaussian_cube import *
from tqdm import tqdm
import pandas as pd 
import numpy as np 
import math
import glob
import os
from scipy.spatial.distance import cdist, squareform, pdist

def get_overlap(R, r, d):
    numerator = np.pi * np.power((R + r - d),2) * (np.power(d,2) + 2*d*r - 3*np.power(r,2) + 2*d*R + 6*r*R - 3*np.power(R,2))
    denomenator = 12*d
    return np.divide(numerator, denomenator) 
def get_vol(r):
    return np.divide(4,3)*np.pi*np.power(r,3)

vdw_radii = {'H':1.2, 'C':1.70, 'N':1.55, 'S':1.8, 'O':1.52,
        'Br':1.85, 'Kr':2.02, 'Cl':1.75, 'Ne':1.54, 'F':1.47}

for mol_func in glob.glob(os.getcwd() + '/ethane/*.cube'):
    mol = Molecule(mol_func)
    crds = mol.coords_and_atoms

    # generate distance matrix for the molecule and save as a pandas dataframe
    crds_array = crds[['x', 'y', 'z']].to_numpy()
    dist_mat = cdist(crds_array, crds_array)
    dist_df = pd.DataFrame(dist_mat, columns = crds['atom'], index = crds['atom'])
    print(dist_df)
    print(dist_df.columns)
    dist_row, dist_col = dist_mat.shape
    overlap_matrix = np.zeros(shape=(dist_row, dist_col))

    idx_cnt = 0

    for index, row in dist_df.iterrows():
        col_cnt = 0
        
        R = vdw_radii[index]

        for el in row:
            
            if el > 0:
                r = vdw_radii[dist_df.columns[col_cnt]]
                if (r+R) < el:
                    overlap_matrix[idx_cnt, col_cnt] = get_overlap(R, r, el)        
            
            col_cnt += 1
        idx_cnt += 1
    tot_vol = 0

    for atom in dist_df.columns:
        tot_vol += get_vol(vdw_radii[atom])

    tot_overlap = np.sum(overlap_matrix)/2
    print(tot_vol)
    print(tot_vol - tot_overlap)
    # print(overlap_matrix)
