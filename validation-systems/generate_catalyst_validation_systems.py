import pandas as pd
import os
import argparse
import smores
from itertools import product
import pathlib
import stk
import rdkit
import re
import logging


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("substituent_csv")
    parser.add_argument("calculation_directory")
    parser.add_argument("catalysts_csv")
    parser.add_argument("-rs", "--replacement_groups", nargs='*', action='store', type=str)
    return parser.parse_args()


def build_stk_molecule(
        parent_group: str,
        substituent: str,
        output_directory: pathlib.Path,
) -> stk.ConstructedMolecule:

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


def gen_molecule_smiles(
        replacement_groups: list,
        parent_group: str,
        substituent: str,
        output_directory: pathlib.Path,
        ) -> str:

    replacement_group_placement_count = 0
    for replacement_group in replacement_groups:
        replacement_group_placement_count += len(re.findall(replacement_group, parent_group))

    try:
        parent_group_smiles = build_stk_molecule(parent_group, substituent, output_directory)

        if replacement_group_placement_count > 1:
            parent_group_smiles = parent_group_smiles.replace("Lu", "Br")
            parent_group_smiles = build_stk_molecule(parent_group_smiles, substituent, output_directory)
    except:
        logging.warning(f'{output_directory} is unparseable by rdkit')
        parent_group_smiles = None

    return parent_group_smiles


def optimize_and_gen_esp(
        substituent_name: str,
        parent_group: str,
        substituent: str,
        replacement_group: list,
        output_directory: pathlib.Path,
) -> None:

    smiles = gen_molecule_smiles(
            replacement_group,
            parent_group,
            substituent,
            output_directory,
            )

    """
    if output_directory.exists():
        pass
    else:
        mol = smores.Molecule(smiles)
        mol.calculate_electrostatic_potential(output_directory, optimize=False)
    """

def main() -> None:
    cli_args = _get_command_line_arguments()

    logging.basicConfig(
            filename=str(pathlib.Path(cli_args.calculation_directory).joinpath(
                'validation_generation_log.log'
                )),
            level=logging.DEBUG
            )

    substituents = pd.read_csv(cli_args.substituent_csv)
    substituents = substituents.iloc[1:]
    catalysts = pd.read_csv(cli_args.catalysts_csv)

    catalyst_smiles = list(catalysts['catalytic_smiles'])

    for molecule_combo in product(catalyst_smiles, list(substituents['substituent_smiles'])):
        substituent_name = substituents.loc[substituents['substituent_smiles'] == molecule_combo[1]]['substituent_name'].values[0]
        catalyst_name = catalysts.loc[catalysts['catalytic_smiles'] == molecule_combo[0]]['catalyst_name'].values[0]
        reaction_name = catalysts.loc[catalysts['catalytic_smiles'] == molecule_combo[0]]['catalytic_reaction_name'].values[0]

        output_directory = pathlib.Path(
                f'{cli_args.calculation_directory}/{reaction_name}/{catalyst_name}/{substituent_name}/'
                )
        output_directory.mkdir(exist_ok=True, parents=True) 
        optimize_and_gen_esp(
                substituent_name,
                molecule_combo[0],
                molecule_combo[1],
                cli_args.replacement_groups,
                output_directory,
            )
        # print('---------------------------------')


if __name__ == '__main__':
    main()
