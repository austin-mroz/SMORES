#!python
import argparse
import json
import pathlib
import sqlite3
import sys
import typing
from dataclasses import dataclass

import atomlite
import rdkit.Chem as rdkit

import smores


def main() -> None:
    args = _get_command_line_arguments()
    args.output_directory.mkdir(parents=True, exist_ok=True)

    database = atomlite.Database(args.database)

    for system in database.get_entries():
        print(system.key)
        entry = atomlite.PropertyEntry(
            key=system.key,
            properties={
                "xyz_file": _replace_directory(
                    pathlib.Path(system.properties["xyz_file"]),
                    args.output_directory,
                ),
                "esp_file": _replace_directory(
                    pathlib.Path(system.properties["esp_file"]),
                    args.output_directory,
                ),
            },
        )
        database.update_properties(entry)
    database.connection.commit()


def _replace_directory(
    input_directory: pathlib.Path,
    output_directory: pathlib.Path,
) -> str:
    input_parts = list(input_directory.parts)[-3:]
    if "xyz" in input_parts[-1]:
        return str(output_directory / input_parts[1] / input_parts[2])
    elif "cube" in input_parts[-1]:
        return str(
            output_directory / input_parts[0] / input_parts[1] / input_parts[2]
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
        default=pathlib.Path.cwd() / "smores_workflow_systems.db",
    )

    parser.add_argument(
        "-o",
        "--output_directory",
        help="The directory into which the results are written.",
        type=pathlib.Path,
        default=_get_output_directory(),
    )
    return parser.parse_args()


def _get_output_directory() -> pathlib.Path:
    return pathlib.Path(str(pathlib.Path.cwd()).replace("work", "data"))


if __name__ == "__main__":
    main()
