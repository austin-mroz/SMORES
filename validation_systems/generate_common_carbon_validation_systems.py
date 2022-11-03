import pandas as pd
import os
import argparse
import smores
from itertools import product
import pathlib

def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("substituent_csv")
    parser.add_argument("calculation_directory")
    parser.add_argument("-r","--replacement_group",type=str)
    parser.set_defaults(replacement_group='Br')
    parser.add_argument("-r2","--replacement_group2",type=str)
    parser.set_defaults(replacement_group='Lu')
    return parser.parse_args()


def gen_molecule_smiles(replacement_group: str,
        r_group: str,
        substituent: str) -> str:
    return substituent.replace(replacement_group, r_group)


def optimize_and_gen_esp(substituent_name: str,
        r_group: str,
        substituent: str,
        replacement_group: str,
        calculation_directory: str) -> None:

    mol_smiles = gen_molecule_smiles(
            replacement_group,
            r_group,
            substituent,
            )
    print(mol_smiles)
    output_directory = pathlib.Path(f'{calculation_directory}/{r_group}/{substituent_name}/')
    if output_directory.exists():
        pass
    else:
        mol = smores.Molecule(mol_smiles)
        mol.calculate_electrostatic_potential(output_directory, optimize=True)


def main() -> None:
    cli_args = _get_command_line_arguments()

    substituents = pd.read_csv(cli_args.substituent_csv)
    
    r_groups = pd.DataFrame([['C', 'c1ccccc1']])
    
    for molecule_combo in product(r_groups, list(substituents['substituent_smiles'])):
        substituent_name = substituents.loc[substituents['substituent_smiles'] == molecule_combo[1]]['substituent_name'].values[0]
        optimize_and_gen_esp(
                substituent_name,
                molecule_combo[0],
                molecule_combo[1],
                cli_args.replacement_group,
                cli_args.calculation_directory,
            )
        print('---------------------------------')


if __name__ == '__main__':
    main()

