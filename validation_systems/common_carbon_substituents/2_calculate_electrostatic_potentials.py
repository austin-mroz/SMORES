import argparse
import csv
import pathlib
import typing
from dataclasses import dataclass

import smores


def main() -> None:
    args = _get_command_line_arguments()
    args.output_directory.mkdir(parents=True, exist_ok=True)
    molecules = tuple(_read_molecules(args.input_file))
    for molecule_input in molecules:
        calculation_directory = args.output_directory / molecule_input.name
        if not calculation_directory.exists():
            molecule = smores.Molecule(molecule_input.smiles)
            molecule.calculate_electrostatic_potential(
                calculation_directory,
                optimize=True,
            )


@dataclass(frozen=True, slots=True)
class Molecule:
    name: str
    smiles: str


def _read_molecules(path: pathlib.Path) -> typing.Iterator[Molecule]:
    with open(path) as csv_file:
        reader = csv.DictReader(csv_file)
        for row in reader:
            yield Molecule(row["name"], row["smiles"])


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Calculate molecular electrostatic potentials.",
    )
    parser.add_argument(
        "-i",
        "--input_file",
        help='A csv file with two columns: "name" and "smiles".',
        type=pathlib.Path,
        default=pathlib.Path(__file__).parent / "1_output" / "smiles.csv",
    )
    parser.add_argument(
        "-o",
        "--output_directory",
        help="The directory into which the results are written.",
        type=pathlib.Path,
        default=pathlib.Path(__file__).parent / "2_output",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
