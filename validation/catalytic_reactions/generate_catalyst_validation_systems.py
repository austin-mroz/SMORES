#!python

import argparse
import logging
import pathlib
import re
from itertools import product

import pandas as pd
import rdkit
import stk

import smores


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("substituent_csv")
    parser.add_argument("calculation_directory")
    parser.add_argument("catalysts_csv")
    parser.add_argument(
        "-rs", "--replacement_groups", nargs="*", action="store", type=str
    )
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
            repeating_unit="AB",
            num_repeating_units=1,
            optimizer=stk.Collapser(scale_steps=False),
        ),
    )
    mol_file_name = output_directory.joinpath("stk_constructed_molecule.mol")

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
        replacement_group_placement_count += len(
            re.findall(replacement_group, parent_group)
        )

    try:
        parent_group_smiles = build_stk_molecule(
            parent_group, substituent, output_directory
        )

        if replacement_group_placement_count > 1:
            parent_group_smiles = parent_group_smiles.replace("Lu", "Br")
            parent_group_smiles = build_stk_molecule(
                parent_group_smiles, substituent, output_directory
            )
    except:
        logging.warning(f"{output_directory} is unparseable by rdkit")
        parent_group_smiles = None

    return parent_group_smiles


def generate_xtb_xyz(
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

    if smiles is not None:
        mol = smores.Molecule(smiles)
        mol.generate_xtb_starting_structure(output_directory)
    else:
        logging.error(f"{output_directory} xyz not generated.")
    # mol.calculate_electrostatic_potential(output_directory, optimize=False)


def generate_xtb_submission(
    calculation_directory: pathlib.Path,
    initial_xyz_file_name: str = "initial_structure.xyz",
    calculation_type: str = "opt",  # options are esp or opt
) -> None:
    xtb_sub = []
    xtb_sub.append("#!/bin/bash")
    xtb_sub.append("")
    xtb_sub.append("run_xtb () {")
    xtb_sub.append("")
    xtb_sub.append("    cd $2")
    xtb_sub.append('    DIR="$(dirname "$1")"')
    xtb_sub.append('    cd "$2$DIR"')
    xtb_sub.append(f'    xyz="{initial_xyz_file_name}"')
    xtb_sub.append(f"    xtb $xyz --{calculation_type}")
    xtb_sub.append('    cd "$2"')
    xtb_sub.append("}")
    xtb_sub.append("")
    xtb_sub.append("export -f run_xtb")
    xtb_sub.append(f'parentDIR="{str(calculation_directory)}"')
    xtb_sub.append("cd $parentDIR")

    xtb_sub_content = "\n".join(xtb_sub)
    with open(f"{calculation_type}_with_xtb.sh", "w") as xtb_sub_file:
        xtb_sub_file.write(f"{xtb_sub_content}\n")
        xtb_sub_file.write(
            f'find . -name "{initial_xyz_file_name}" -exec bash -c '
        )
        xtb_sub_file.write("'")
        xtb_sub_file.write('run_xtb "$@" "')
        xtb_sub_file.write("'${parentDIR}'")
        xtb_sub_file.write('"')
        xtb_sub_file.write("'")
        xtb_sub_file.write(" bash {} \;")


def main() -> None:
    cli_args = _get_command_line_arguments()

    logging.basicConfig(
        filename="validation_generation_log.log",
        level=logging.DEBUG,
    )

    substituents = pd.read_csv(cli_args.substituent_csv)
    substituents = substituents.iloc[1:]
    catalysts = pd.read_csv(cli_args.catalysts_csv)

    catalyst_smiles = list(catalysts["catalytic_smiles"])

    for molecule_combo in product(
        catalyst_smiles, list(substituents["substituent_smiles"])
    ):
        substituent_name = substituents.loc[
            substituents["substituent_smiles"] == molecule_combo[1]
        ]["substituent_name"].values[0]
        catalyst_name = catalysts.loc[
            catalysts["catalytic_smiles"] == molecule_combo[0]
        ]["catalyst_name"].values[0]
        reaction_name = catalysts.loc[
            catalysts["catalytic_smiles"] == molecule_combo[0]
        ]["catalytic_reaction_name"].values[0]

        output_directory = pathlib.Path(
            f"{cli_args.calculation_directory}/{reaction_name}/{catalyst_name}/{substituent_name}/"
        )
        print(output_directory)

        output_directory.mkdir(exist_ok=True, parents=True)

        generate_xtb_xyz(
            substituent_name,
            molecule_combo[0],
            molecule_combo[1],
            cli_args.replacement_groups,
            output_directory,
        )
    generate_xtb_submission(
        cli_args.calculation_directory,
    )
    # print('---------------------------------')


if __name__ == "__main__":
    main()
