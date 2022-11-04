import argparse
import csv
import pathlib
import sys
import typing
from dataclasses import dataclass

import smores


def main() -> None:
    args = _get_command_line_arguments()
    args.output_directory.mkdir(parents=True, exist_ok=True)
    molecules = tuple(_read_molecules(args.input_file))
    csv_output = args.output_directory / "xyz_files.csv"
    for molecule_input in molecules:
        calculation_directory = args.output_directory / molecule_input.name
        if not calculation_directory.exists():

            try:
                molecule = smores.Molecule(molecule_input.smiles)
                molecule.calculate_electrostatic_potential(
                    calculation_directory,
                    optimize=True,
                )
                _append_to_csv(
                    path=csv_output,
                    xyz_file=calculation_directory / "geom.xyz",
                    dummy_index=molecule_input.dummy_index,
                    attached_index=molecule_input.attached_index,
                )

            except Exception as ex:
                print(
                    f"Issue with: {molecule_input.name}\n{ex}",
                    file=sys.stderr,
                )


@dataclass(frozen=True, slots=True)
class Molecule:
    name: str
    smiles: str
    dummy_index: int
    attached_index: int


def _read_molecules(path: pathlib.Path) -> typing.Iterator[Molecule]:
    with open(path) as csv_file:
        reader = csv.DictReader(csv_file)
        for row in reader:
            yield Molecule(
                name=row["name"],
                smiles=row["smiles"],
                dummy_index=int(row["dummy_index"]),
                attached_index=int(row["attached_index"]),
            )


def _append_to_csv(
    path: pathlib.Path,
    xyz_file: pathlib.Path,
    dummy_index: int,
    attached_index: int,
) -> None:

    with open(path, "a") as csv_file:
        writer = csv.DictWriter(
            csv_file,
            fieldnames=["xyz_file", "dummy_index", "attached_index"],
        )
        writer.writeheader()
        writer.writerow(
            {
                "xyz_file": str(xyz_file.resolve()),
                "dummy_index": str(dummy_index),
                "attached_index": str(attached_index),
            }
        )


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
