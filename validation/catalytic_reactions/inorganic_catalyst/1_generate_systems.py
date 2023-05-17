#!python
import argparse
import csv
import pathlib

import atomlite
import rdkit.Chem as rdkit

import smores


def main() -> None:
    args = _get_command_line_arguments()
    args.output_directory.mkdir(parents=True, exist_ok=True)

    cores = {
        smores.rdkit_from_smiles(
            "CC1=[N](C2=C(C3C4=CC=CC=C4C5C6=CC=CC=C63)C5=C(OC)C7=C2[C@@H]8C9=CC=CC=C9[C@H]7C%10=CC=CC=C%108)[Ni]([Lu])([Lu])[N](Br)=C1C"
        ): "Ni_catalyst",
    }

    # naming scheme is R1_R1

    substituents = {
        smores.rdkit_from_smiles(
            "BrC1=C(C(C2=CC=C(C)C=C2)C3=CC=C(C)C=C3)C=C(C)C=C1C(C4=CC=C(C)C=C4)C5=CC=C(C)C=C5"
        ): "Me_Me",
        smores.rdkit_from_smiles(
            "BrC1=C(C(C2=CC=C(F)C=C2)C3=CC=C(F)C=C3)C=C(C)C=C1C(C4=CC=C(F)C=C4)C5=CC=C(F)C=C5"
        ): "Me_F",
        smores.rdkit_from_smiles(
            "BrC1=C(C(C2=CC=C([H])C=C2)C3=CC=C([H])C=C3)C=C(C)C=C1C(C4=CC=C([H])C=C4)C5=CC=C([H])C=C5"
        ): "Me_H",
        smores.rdkit_from_smiles(
            "BrC1=C(C(C2=CC=C([H])C=C2)C3=CC=C([H])C=C3)C=C(Cl)C=C1C(C4=CC=C([H])C=C4)C5=CC=C([H])C=C5"
        ): "Cl_H",
        smores.rdkit_from_smiles(
            "BrC1=C(C(C2=CC=C([H])C=C2)C3=CC=C([H])C=C3)C=C(OC)C=C1C(C4=CC=C([H])C=C4)C5=CC=C([H])C=C5"
        ): "MeO_H",
    }

    system_database = atomlite.Database("nickel_catalyst_systems.db")

    database_entries = []

    for combo in smores.combine(cores, substituents, join_atom="Br"):
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
