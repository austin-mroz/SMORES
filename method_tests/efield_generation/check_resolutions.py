import sys; sys.path.insert(0,'/home/he/bin/github_clones/STREUSEL/streusel/')
from streusel import gaussian_cube

from molecule import *

import matplotlib.pyplot as plt
from tqdm import tqdm
import pandas as pd
import glob
import os

vols = pd.DataFrame()

for res in tqdm(range(50, 200, 20)):
    calc_path = f'{os.getcwd()}/res_check/res_{res}/'
    mol = Molecule(os.getcwd(), calc_path)
    mol.init_from_xyz(f'{os.getcwd()}/ethane_opt.xyz')

    mol.gen_esp_cube(res=res)

    esp = f'{os.getcwd()}/res_check/res_{res}/ESP.cube'
    mol = gaussian_cube.Molecule(esp)
    mol.get_efield()
    mol.sample_efield()

    nrow = pd.DataFrame([[res, mol.vol]])
    tempdf = pd.concat([vols, nrow])
    vols = pd.DataFrame(tempdf)

print(vols)
vols.to_csv('res_check.csv')

