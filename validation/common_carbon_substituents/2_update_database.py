#!python
import argparse
import json
import pathlib
import sqlite3
import sys
import typing
from dataclasses import dataclass

import rdkit.Chem as rdkit

import smores

# import smores.psi4


def main() -> None:
    args = _get_command_line_arguments()
    args.output_directory.mkdir(parents=True, exist_ok=True)
    database = sqlite3.connect(args.database)
    cursor = database.cursor()
    molecules = tuple(_get_molecules(database))

    for molecule_input in molecules:
        calculation_directory = args.output_directory / molecule_input.name

        cursor.execute(
            "UPDATE molecules SET properties=json_insert(properties,?,?) WHERE key=?",
            (
                "$.esp_file",
                str(calculation_directory / "ESP.cube"),
                molecule_input.name,
            ),
        )

        database.commit()

    database.close()


@dataclass(frozen=True, slots=True)
class Molecule:
    name: str
    xyz_file: pathlib.Path


def _get_molecules(database: sqlite3.Connection) -> typing.Iterator[Molecule]:
    for key, _, properties in database.execute("SELECT * FROM molecules"):
        yield Molecule(
            name=key,
            xyz_file=json.loads(properties)["xyz_file"],
        )


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Calculate molecular electrostatic potentials.",
    )
    parser.add_argument(
        "-d",
        "--database",
        help=(
            'An atomlite database file with properties: "core", "substituent", '
            '"dummy_index", "attached_index", "xyz_file".'
        ),
        type=pathlib.Path,
        default=pathlib.Path.cwd() / "common_carbon_substituents.db",
    )
    parser.add_argument(
        "-o",
        "--output_directory",
        help="The directory into which the results are written.",
        type=pathlib.Path,
        default=pathlib.Path.cwd() / "2_output",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
