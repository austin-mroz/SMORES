#!python
import argparse
import csv
import pathlib

import atomlite
import rdkit.Chem as rdkit
import vanila

import smores


def main() -> None:
    args = _get_command_line_arguments()
    args.output_directory.mkdir(parents=True, exist_ok=True)

    cores = {
        smores.rdkit_from_smiles("BrN1C=CN=C1C2N(Br)CC[NH+]2Br"): "JOC_2",
    }

    # naming scheme is R1_R1

    substituents = {
        smores.rdkit_from_smiles("CBr"): "CH3",
        smores.rdkit_from_smiles("CCBr"): "CH2CH3",
        smores.rdkit_from_smiles("CCCBr"): "CH2CH2CH3",
        smores.rdkit_from_smiles("CCCCBr"): "CH2CH2CH2CH3",
        smores.rdkit_from_smiles("FC(F)(F)CBr"): "CH2CF3",
        smores.rdkit_from_smiles("FC(F)(F)CCBr"): "CH2CH2CF3",
    }

    system_database = atomlite.Database("catalyst_systems.db")

    database_entries = []

    for combo in vanila.combine._get_one_core_n_cap_ions(
        substituents,
        cores,
        join_atom="Br",
        output_directory=pathlib.Path.cwd(),
    ):
        print(combo)
        exit()
        name = f"{cores[combo.core]}_{substituents[combo.substituent]}"
        xyz_file = args.output_directory / f"{name}.xyz"
        rdkit.MolToXYZFile(combo.product, str(xyz_file))

        database_entries.append(
            atomlite.Entry.from_rdkit(
                key=name,
                molecule=combo.product,
                properties={
                    "core": cores[combo.core],
                    "substituent": substituents[combo.substituent],
                    "dummy_index": combo.dummy_index,
                    "attached_index": combo.attached_index,
                    "xyz_file": str(xyz_file),
                },
            )
        )

    system_database.add_entries(database_entries)

    system_database.connection.commit()


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate the structures of systems which are to be validated."
        ),
    )
    parser.add_argument(
        "-o",
        "--output_directory",
        default=pathlib.Path.cwd() / "1_output",
        type=pathlib.Path,
        help="The directory into which the output files are written.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
