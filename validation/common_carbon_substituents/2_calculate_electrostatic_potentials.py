#!python
import argparse
import csv
import pathlib
import sys
import typing
from dataclasses import dataclass

import rdkit.Chem as rdkit

import smores
import smores.psi4

_OUTPUT_CSV_COLUMNS = (
    "name",
    "core",
    "substituent",
    "smiles",
    "xyz_file",
    "dummy_index",
    "attached_index",
)


def main() -> None:
    args = _get_command_line_arguments()
    args.output_directory.mkdir(parents=True, exist_ok=True)
    molecules = tuple(_read_molecules(args.input_file))

    csv_output = args.output_directory / "xyz_files.csv"
    if not csv_output.exists():
        _write_csv_header(csv_output)

    for molecule_input in molecules:
        calculation_directory = args.output_directory / molecule_input.name
        if not calculation_directory.exists():
            try:
                molecule = rdkit.MolFromXYZFile(str(molecule_input.xyz_file))
                smores.psi4.calculate_electrostatic_potential(
                    molecule=molecule,
                    output_directory=calculation_directory,
                    grid_origin=(-5, -5, -5),
                    grid_length=10.0,
                    num_voxels_per_dimension=50,
                    optimize=True,
                    num_threads=16,
                )
                _append_to_csv(
                    path=csv_output,
                    name=molecule_input.name,
                    core=molecule_input.core,
                    substituent=molecule_input.substituent,
                    smiles=molecule_input.smiles,
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
    core: str
    substituent: str
    smiles: str
    dummy_index: int
    attached_index: int
    xyz_file: pathlib.Path


def _read_molecules(path: pathlib.Path) -> typing.Iterator[Molecule]:
    with open(path) as csv_file:
        reader = csv.DictReader(csv_file)
        for row in reader:
            yield Molecule(
                name=row["name"],
                core=row["core"],
                substituent=row["substituent"],
                smiles=row["smiles"],
                dummy_index=int(row["dummy_index"]),
                attached_index=int(row["attached_index"]),
                xyz_file=pathlib.Path(row["xyz_file"]),
            )


def _write_csv_header(path: pathlib.Path) -> None:
    with open(path, "w") as csv_file:
        writer = csv.DictWriter(
            csv_file,
            fieldnames=_OUTPUT_CSV_COLUMNS,
        )
        writer.writeheader()


def _append_to_csv(
    path: pathlib.Path,
    name: str,
    core: str,
    substituent: str,
    smiles: str,
    xyz_file: pathlib.Path,
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
                "name": name,
                "core": core,
                "substituent": substituent,
                "smiles": smiles,
                "xyz_file": xyz_file.resolve(),
                "dummy_index": dummy_index,
                "attached_index": attached_index,
            },
        )


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Calculate molecular electrostatic potentials.",
    )
    parser.add_argument(
        "-i",
        "--input_file",
        help=(
            'A csv file with columns: "name", "core", "substituent", '
            '"smiles", "dummy_index", "attached_index", "xyz_file".'
        ),
        type=pathlib.Path,
        default=pathlib.Path.cwd() / "1_output" / "xyz_files.csv",
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
