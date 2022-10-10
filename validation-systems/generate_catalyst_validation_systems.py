import pandas as pd
import os
import argparse
import smores
from itertools import product
import pathlib
import stk
import rdkit

def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("substituent_csv")
    parser.add_argument("calculation_directory")
    parser.add_argument("catalysts_csv") 
    parser.add_argument("-rs","--replacement_groups",nargs='*',action='store',type=str)
    return parser.parse_args()


def gen_molecule_smiles(
        replacement_groups: str,
        parent_group: str,
        substituent: str,
        output_directory: pathlib.Path,
        ) -> str:

    parent_mol = stk.BuildingBlock(
            smiles=parent_group,
            functional_groups=[stk.BromoFactory()],
            )
    substituent_mol = stk.BuildingBlock(
            smiles=substituent,
            functional_groups=[stk.BromoFactory()],
            )
    full_catalyst = stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(parent_mol, substituent_mol),
                repeating_unit='AB',
                num_repeating_units=1,
                optimizer=stk.Collapser(scale_steps=False),
                ),
            )
    mol_file_name = output_directory.joinpath('stk_constructed_molecule.mol')

    full_catalyst.write(str(mol_file_name))

    rdmol = rdkit.Chem.rdmolfiles.MolFromMolFile(str(mol_file_name))

    return rdkit.Chem.rdmolfiles.MolToSmiles(rdmol)


def optimize_and_gen_esp(
        substituent_name: str,
        parent_group: str,
        substituent: str,
        replacement_group: str,
        calculation_directory: str) -> None:

    smiles = gen_molecule_smiles(
            replacement_group,
            parent_group,
            substituent,
            pathlib.Path(calculation_directory),
            )

    output_directory = pathlib.Path(f'{calculation_directory}/{parent_group}/{substituent_name}/')
    if output_directory.exists():
        pass
    else:
        mol = smores.Molecule(smiles)
        mol.calculate_electrostatic_potential(output_directory, optimize=True)
        exit()

def main() -> None:
    cli_args = _get_command_line_arguments()

    substituents = pd.read_csv(cli_args.substituent_csv)
    substituents = substituents.iloc[1:]
    catalysts = pd.read_csv(cli_args.catalysts_csv)

    catalyst_smiles = list(catalysts['catalytic_smiles'])

    for molecule_combo in product(catalyst_smiles, list(substituents['substituent_smiles'])):
        substituent_name = substituents.loc[substituents['substituent_smiles'] == molecule_combo[1]]['substituent_name'].values[0]
        print(substituent_name)
        print(molecule_combo[1])
        print(molecule_combo[0])
        optimize_and_gen_esp(
                substituent_name,
                molecule_combo[0],
                molecule_combo[1],
                cli_args.replacement_groups,
                cli_args.calculation_directory,
            )
        print('---------------------------------')


if __name__ == '__main__':
    main()

