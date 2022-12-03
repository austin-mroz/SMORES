#!python
import argparse
import csv
import pathlib
import typing
from dataclasses import dataclass

import statsmodels.api as sm

_OUTPUT_CSV_COLUMNS = (
    "name",
    "core",
    "substituent",
    "smiles",
    "xyz_file",
    "dummy_index",
    "attached_index",
    "radii_type",
    "L",
    "B1",
    "B5",
)


@dataclass(frozen=True, slots=True)
class SterimolFit:
    name: str
    core: str
    substituent: str
    radii_type: str
    L: float
    B1: float
    B5: float
    experimental_ddG: float
    predicted_ddG: float


@dataclass(frozen=True, slots=True)
class CsvRow:
    name: str
    core: str
    substituent: str
    smiles: str
    xyz_file: pathlib.Path
    dummy_index: int
    attached_index: int
    radii_type: str
    L: float
    B1: float
    B5: float


def main() -> None:
    args = _get_command_line_arguments()
    args.output_directory.mkdir(parents=True, exist_ok=True)

    sterimol_parameters = tuple(_get_rows(args.sterimol_parameters))
    experimental_results = tuple(_get_rows(args.experimental_results))


def _fit_sterimol_parameters(
    sterimol_parameters: typing.Iterator[CsvRow],
    experimental_results: typing.Iterator[CsvRow],
    catalysis_reaction: str,
) -> typing.Iterator[SterimolFit]:
    pass


def _get_rows(
    path: pathlib.Path,
) -> typing.Iterator[CsvRow]:

    with open(path) as csv_file:
        reader = csv.DictReader(csv_file)
        for row in reader:
            yield CsvRow(
                name=row["name"],
                core=row["core"],
                substituent=row["substituent"],
                smiles=row["smiles"],
                xyz_file=pathlib.Path(row["xyz_file"]),
                dummy_index=int(row["dummy_index"]),
                attached_index=int(row["attached_index"]),
                radii_type=row["radii_type"],
                L=float(row["L"]),
                B1=float(row["B1"]),
                B5=float(row["B5"]),
            )


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fit sterimol parameters to experimental results",
    )

    parser.add_argument(
        "-s",
        "--sterimol_parameters",
        help=("A csv file containing sterimol parameters"),
        type=pathlib.Path,
        default=pathlib.Path.cwd()
        .joinpath("4_output")
        .joinpath("steric_parameters.csv"),
    )
    parser.add_argument(
        "-e",
        "--experimental_results",
        help=("a csv file containing experimental ddG"),
        type=pathlib.Path,
    )
    parser.add_argument(
        "-o",
        "--output_directory",
        type=pathlib.Path,
        default=pathlib.Path.cwd() / "6_output",
    )

    parser.add_argument(
        "-c",
        "--catalysis_reaction",
        help=("Catalysis reaction to fit"),
        type=str,
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--from_radii", action="store_true")
    group.add_argument("--from_cube", action="store_true")

    return parser.parse_args()


if __name__ == "__main__":
    main()
