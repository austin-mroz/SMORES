import argparse
import logging
import pathlib
import smores
import streusel
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

        r_group_name = pathlib.PurePath(xyz_path).parts[-3]
        substituent_name = pathlib.PurePath(xyz_path).parts[-2]

        substituent_smiles = substituent_dict.loc[substituent_dict['substituent_name'] == substituent_name]['substituent_smiles'].values[0]

        mol_path = convert_xyz_to_mol(xyz_path)

        # use rdkit fragment to get the part of the molecule that we care about
        full_mol = rdkit.MolFromMolFile(str(mol_path), removeHs=False)

        full_mol_list.append(full_mol)
        
    core = rdkit.MolFromSmiles(r_group_name)

    full_mols = full_mol_list # deepcopy(full_mol_list)
    pattern = core
    
    # find matching atoms in each molecule -- this should be the core
    matching = [molecule.GetSubstructMatch(pattern) for molecule in full_mols]

    img2 = Draw.MolsToGridImage(
            full_mols,
            molsPerRow=3,
            )
    img2.save('grid2.png')

    img = Draw.MolsToGridImage(
            full_mols,
            molsPerRow=3,
            highlightAtomLists=matching,
            )
    print(matching)

    img.save('grid.png')


if __name__ == '__main__':
    main()
