import argparse
import logging
import pathlib
import smores
import morfeus
import dbstep.Dbstep as db
import seaborn as sns
import pandas as pd
import glob
from openbabel import openbabel
import rdkit.Chem.AllChem as rdkit
from rdkit.Chem import rdDepictor, Draw
from rdkit.Chem import rdRGroupDecomposition
from visualize_core import highlight_rgroups
import numpy as np
from itertools import product


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("calculation_directory")
    parser.add_argument("substituent_dictionary_csv")
    parser.add_argument("output_directory")
    return parser.parse_args()


def convert_xyz_to_mol(
        xyz_path: pathlib.Path | str,
) -> pathlib.Path:

    if isinstance(xyz_path, str):
        xyz_path = pathlib.Path(xyz_path)
    mol_path = xyz_path.parent.joinpath('geom.mol')

    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("xyz", "mol")

    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, str(xyz_path))
    obConversion.WriteFile(mol, str(mol_path))

    return mol_path


def draw_it(full_mol, groups, idxs, qcore):
    return highlight_rgroups(full_mol,groups,qcore,lbls=('R1')) # ,'R2','R3','R4'))


def main() -> None:
    cli_args = _get_command_line_arguments()

    total_morfeus_sterimol_df = pd.DataFrame()

    substituent_dict = pd.read_csv(cli_args.substituent_dictionary_csv)

    full_mol_list = []

    for calculation_directory in pathlib.Path(cli_args.calculation_directory).rglob('geom.xyz'):
        xyz_path = pathlib.Path(calculation_directory)
        print(xyz_path)

        r_group_smiles = pathlib.PurePath(xyz_path).parts[-3]
        substituent_name = pathlib.PurePath(xyz_path).parts[-2]

        substituent_smiles = substituent_dict.loc[substituent_dict['substituent_name'] == substituent_name]['substituent_smiles'].values[0]

        mol_path = convert_xyz_to_mol(xyz_path)

        # use rdkit fragment to get the part of the molecule that we care about
        full_mol = rdkit.MolFromMolFile(str(mol_path), removeHs=False)

        core_pattern = rdkit.MolFromSmiles(r_group_smiles)
        subt_pattern = rdkit.AddHs(rdkit.MolFromSmiles(substituent_smiles.replace('Br',"")))

        core_match = full_mol.GetSubstructMatch(core_pattern)
        subt_match = full_mol.GetSubstructMatch(subt_pattern)

        core_coords = []
        core_elements = []
        # core coordinates
        for s in core_match:
            core_elements.append(full_mol.GetAtoms()[s].GetSymbol())
            x, y, z = list(full_mol.GetConformer().GetAtomPosition(s))
            core_coords.append([x, y, z])
        core_xyz = smores.utilities.XyzData(
                elements=core_elements,
                coordinates=core_coords,
                )
        
        smores.utilities.write_xyz('core_xyz.xyz', core_xyz)
        
        full_mol_xyz = smores.utilities.read_xyz(calculation_directory)
        
        core_xyz_array = smores.utilities.get_full_xyz_array(core_xyz)
        full_mol_xyz_array = smores.utilities.get_full_xyz_array(
                full_mol_xyz,
                removeHs=True,
                )

        full_mol_C_coords = np.round(np.delete(full_mol_xyz_array, 0, axis=1).astype(np.float32),3)
        core_xyz_C_coords = np.round(np.delete(core_xyz_array, 0, axis=1).astype(np.float32),3)

        full_C_coords_df = pd.DataFrame(full_mol_C_coords)
        core_C_coords_df = pd.DataFrame(core_xyz_C_coords)

        subst_coords_df = pd.concat([full_C_coords_df, core_C_coords_df]).drop_duplicates(keep=False)
        subst_coords = np.asarray(subst_coords_df)

        for core_coord, subt_coord in product(core_xyz_C_coords, subst_coords):
            if np.linalg.norm(core_coord.reshape(-1,1).T-subt_coord) < 1.6:
                core_min_C_coord = core_coord
                subt_min_C_coord = subt_coord
                break
        print(core_min_C_coord)
        print(subt_min_C_coord)
        print(core_min_C_coord == np.round(full_mol_xyz.coordinates,3))
        core_connecting_atom_index = np.where(np.all(core_min_C_coord == np.round(full_mol_xyz.coordinates,3),axis=1))
        sybt_connecting_atom_index = np.where(np.all(subt_min_C_coord == np.round(full_mol_xyz.coordinates,3),axis=1))
        exit()


        exit()
        print(np.where(np.all(core_min_C_coord==full_mol_xyz,axis=1)))
        exit()


        print('substituent match')
        for s in subt_match:
            print(full_mol.GetAtoms()[s].GetSymbol())
            print(list(full_mol.GetConformer().GetAtomPositions(s)))
        print(substituent_smiles.replace('Br',""))


if __name__ == '__main__':
    main()
