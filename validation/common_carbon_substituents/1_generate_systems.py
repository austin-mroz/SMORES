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
        smores.rdkit_from_smiles("BrC"): "C",
        smores.rdkit_from_smiles("Brc1ccccc1"): "c1ccccc1",
    }
    substituents = {
        smores.rdkit_from_smiles("Br"): "H",
        smores.rdkit_from_smiles("BrC"): "Me",
        smores.rdkit_from_smiles("BrCC"): "Et",
        smores.rdkit_from_smiles("Brc1ccccc1"): "Ph",
        smores.rdkit_from_smiles("BrCc1ccccc1"): "Bn",
        smores.rdkit_from_smiles("BrC(C)CC"): "CH2-iPr",
        smores.rdkit_from_smiles("BrC(C)(C)CC"): "CH2-tBu",
        smores.rdkit_from_smiles("CCCC(Br)CCC"): "CHPr2",
        smores.rdkit_from_smiles("N#CBr"): "Cy",
        smores.rdkit_from_smiles("CC(C)C(Br)C(C)C"): "CH(i-Pr)2",
        smores.rdkit_from_smiles("CCC(Br)CC"): "CHEt2",
        smores.rdkit_from_smiles("CCC(Br)(CC)CC"): "CEt3",
    }

    system_database = atomlite.Database("common_carbon_substituents.db")

    database_entries = []

    for combo in smores.combine(cores, substituents):
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
