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
import smores.psi4


def main() -> None:
    args = _get_command_line_arguments()
    args.output_directory.mkdir(parents=True, exist_ok=True)
    database = sqlite3.connect(args.database)
    cursor = database.cursor()
    molecules = tuple(_get_molecules(database))

    for molecule_input in molecules:
        calculation_directory = args.output_directory / molecule_input.name
        if not calculation_directory.exists():
            try:
                print(molecule_input.name)
                molecule = rdkit.MolFromXYZFile(str(molecule_input.xyz_file))
                smores.psi4.calculate_electrostatic_potential(
                    molecule=molecule,
                    output_directory=calculation_directory,
                    grid_origin=(-20, -20, -20),
                    grid_length=26.0,
                    num_voxels_per_dimension=100,
                    optimize=True,
                    num_threads=20,
                )

                cursor.execute(
                    "UPDATE molecules SET properties=json_insert(properties,?,?) WHERE key=?",
                    (
                        "$.esp_file",
                        str(calculation_directory / "ESP.cube"),
                        molecule_input.name,
                    ),
                )

                database.commit()
            except Exception as ex:
                print(
                    f"Issue with: {molecule_input.name}\n{ex}",
                    file=sys.stderr,
                )

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
        default=pathlib.Path.cwd() / "catalyst_systems.db",
    )
    parser.add_argument(
        "-o",
        "--output_directory",
        default=_get_output_directory(),
        type=pathlib.Path,
        help="The directory into which the output files are written.",
    )
    return parser.parse_args()


def _get_output_directory() -> pathlib.Path:
    return pathlib.Path(
        str(pathlib.Path.cwd() / "2_output").replace("work", "data")
    )


if __name__ == "__main__":
    main()
