#!python

import argparse
import csv
import pathlib
import sys
import typing
from collections import abc
from dataclasses import dataclass

import rdkit.Chem as rdkit

import smores
import smores.psi4

_OUTPUT_CSV_COLUMNS = (
    "reaction_name",
    "name",
    "core",
    "substituent",
    "smiles",
    "xyz_file",
    "fragments_file",
    "dummy_index",
    "attached_index",
)


def main() -> None:
    args = _get_command_line_arguments()
    args.output_directory.mkdir(parents=True, exist_ok=True)

    for input_file in args.input_file:
        _calculate_electrostatic_potentials(
            molecules=tuple(_read_molecules(input_file)),
            output_directory=args.output_directory,
        )


@dataclass(frozen=True, slots=True)
class Molecule:
    reaction_name: str
    name: str
    core: str
    substituent: str
    smiles: str
    dummy_index: int
    attached_index: int
    xyz_file: pathlib.Path
    fragments_file: pathlib.Path


def _calculate_electrostatic_potentials(
    molecules: abc.Collection[Molecule],
    output_directory: pathlib.Path,
) -> None:

    for molecule in molecules:
        reaction_directory = output_directory / molecule.reaction_name
        reaction_directory.mkdir(parents=True, exist_ok=True)
        csv_path = reaction_directory / "xyz_files.csv"

        with open(csv_path, "w") as csv_file:
            writer = csv.DictWriter(
                csv_file,
                fieldnames=_OUTPUT_CSV_COLUMNS,
            )
            writer.writeheader()

    for molecule_input in molecules:
        reaction_directory = output_directory / molecule_input.reaction_name
        csv_path = reaction_directory / "xyz_files.csv"

        calculation_directory = reaction_directory / molecule_input.name

        if not calculation_directory.exists():
            calculation_directory.mkdir(parents=True)
            try:
                molecule = rdkit.MolFromXYZFile(str(molecule_input.xyz_file))
                smores.psi4.calculate_electrostatic_potential(
                    molecule=molecule,
                    output_directory=calculation_directory,
                    grid_origin=(-5, -5, -5),
                    grid_length=10.0,
                    num_voxels_per_dimension=50,
                    optimize=False,
                    num_threads=20,
                )
                _append_to_csv(
                    path=csv_path,
                    reaction_name=molecule_input.reaction_name,
                    name=molecule_input.name,
                    core=molecule_input.core,
                    substituent=molecule_input.substituent,
                    smiles=molecule_input.smiles,
                    xyz_file=calculation_directory / "xtbopt.xyz",
                    fragments_file=molecule_input.fragments_file,
                    dummy_index=molecule_input.dummy_index,
                    attached_index=molecule_input.attached_index,
                )

            except Exception as ex:
                print(
                    f"Issue with: {molecule_input.name}\n{ex}",
                    file=sys.stderr,
                )


def _append_to_csv(
    path: pathlib.Path,
    reaction_name: str,
    name: str,
    core: str,
    substituent: str,
    smiles: str,
    xyz_file: pathlib.Path,
    fragments_file: pathlib.Path,
    dummy_index: int,
    attached_index: int,
) -> None:
    with open(path, "a") as csv_file:
        writer = csv.DictWriter(
            csv_file,
            fieldnames=_OUTPUT_CSV_COLUMNS,
        )
        writer.writerow(
            {
                "reaction_name": reaction_name,
                "name": name,
                "core": core,
                "substituent": substituent,
                "smiles": smiles,
                "xyz_file": xyz_file.resolve(),
                "fragments_file": fragments_file,
                "dummy_index": dummy_index,
                "attached_index": attached_index,
            },
        )


def _read_molecules(input_file: pathlib.Path) -> typing.Iterator[Molecule]:
    with open(input_file) as csv_file:
        reader = csv.DictReader(csv_file)
        for row in reader:
            yield Molecule(
                reaction_name=row["reaction_name"],
                name=row["name"],
                core=row["core"],
                substituent=row["substituent"],
                smiles=row["smiles"],
                dummy_index=int(row["dummy_index"]),
                attached_index=int(row["attached_index"]),
                xyz_file=pathlib.Path(row["xyz_file"]),
                fragments_file=pathlib.Path(row["fragments_file"]),
            )


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Calculate molecular electrostatic potentials.",
    )
    parser.add_argument(
        "-i",
        "--input_file",
        help=(
            'A csv file with columns: "reaction_name", "name", "core", '
            '"substituent", '
            '"smiles", "dummy_index", "attached_index", "xyz_file".'
        ),
        type=pathlib.Path,
        default=pathlib.Path.cwd()
        .joinpath("2_output")
        .glob("*/xyz_files.csv"),
        nargs="+",
    )
    parser.add_argument(
        "-o",
        "--output_directory",
        help="The directory into which the results are written.",
        type=pathlib.Path,
        default=pathlib.Path.cwd() / "3_output",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
