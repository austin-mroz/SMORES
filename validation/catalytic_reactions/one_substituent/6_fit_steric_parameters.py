#!python
import argparse
import pathlib
import typing
from dataclasses import dataclass

import pandas as pd

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


def main() -> None:
    args = _get_command_line_arguments()
    args.output_directory.mkdir(parents=True, exist_ok=True)

    sterimol_parameters = pd.read_csv(args.sterimol_parameters)
    experimental_results = pd.read_csv(args.experimental_results)
    for core in sterimol_parameters["core"].unique():
        _fit_sterimol_parameters(
            sterimol_parameters=sterimol_parameters,
            experimental_results=experimental_results,
            core=core,
        )


def _fit_sterimol_parameters(
    sterimol_parameters: pd.DataFrame,
    experimental_results: pd.DataFrame,
    core: str,
) -> typing.Iterator[SterimolFit]:

    parameter_combinations = _powerset("L", "B1", "B5")
    for parameter_combination in parameter_combinations:
        experimental_results[]



def _powerset(*items: str) -> list[list[str]]:
    subsets: list[list[str]] = [[]]
    for item in items:
        subsets += [subset + [item] for subset in subsets]
    # Do not return the empty set.
    return subsets[1:]


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fit sterimol parameters to experimental results.",
    )

    parser.add_argument(
        "--sterimol_parameters",
        help="A csv file containing sterimol parameters.",
        type=pathlib.Path,
        default=pathlib.Path.cwd() / "4_output" / "steric_parameters.csv",
    )
    parser.add_argument(
        "--experimental_results",
        help="A csv file containing experimental ddG.",
        type=pathlib.Path,
        default=pathlib.Path.cwd() / "exprimental_ddG.csv"
    )
    parser.add_argument(
        "-o",
        "--output_directory",
        type=pathlib.Path,
        default=pathlib.Path.cwd() / "6_output",
    )

    return parser.parse_args()


if __name__ == "__main__":
    main()
