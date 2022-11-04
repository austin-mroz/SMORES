import argparse
import glob
import logging
import pathlib

import dbstep.Dbstep as db
import morfeus
import numpy
import pandas as pd
import seaborn as sns
import streusel

import smores
from smores.calculators import SterimolParameter


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("calculation_directory")
    parser.add_argument("output_directory")
    return parser.parse_args()


def get_morfeus_sterimol_dataframe(
    sterimol_parameters: list[SterimolParameter],
    substituent_name: str,
    r_group_name: str,
) -> pd.DataFrame:

    morfeus_sterimol_df = pd.DataFrame()

    for sterimol_ in sterimol_parameters[0]:

        nrowLval = pd.DataFrame(
            [
                [
                    substituent_name,
                    r_group_name,
                    sterimol_.radii_type,
                    "L_value",
                    sterimol_.L,
                ]
            ]
        )
        nrowB1val = pd.DataFrame(
            [
                [
                    substituent_name,
                    r_group_name,
                    sterimol_.radii_type,
                    "B1_value",
                    sterimol_.B1,
                ]
            ]
        )
        nrowB5val = pd.DataFrame(
            [
                [
                    substituent_name,
                    r_group_name,
                    sterimol_.radii_type,
                    "B5_value",
                    sterimol_.B5,
                ]
            ]
        )
        tempdf = pd.concat(
            [morfeus_sterimol_df, nrowLval, nrowB1val, nrowB5val]
        )
        morfeus_sterimol_df = pd.DataFrame(tempdf)
    morfeus_sterimol_df.columns = [
        "substituent",
        "r_group",
        "radii_type",
        "sterimol parameter",
        "value",
    ]

    return morfeus_sterimol_df


def write_sterimol_to_output(
    sterimol_parameters: list[SterimolParameter],
    output_path: pathlib.Path,
    r_group_name: str,
    substituent_name: str,
) -> None:
    output_lines = []
    output_lines.append("method, radii_type, L, B1, B5")
    for sterimol_ in sterimol_parameters[0]:
        print(sterimol_)
        output_lines.append(
            f"{sterimol_.method}, {sterimol_.radii_type}, {sterimol_.L}, {sterimol_.B1}, {sterimol_.B5}"
        )
    content = "\n".join(output_lines)
    with open(
        output_path.joinpath(f"{r_group_name}_{substituent_name}.csv"),
        "w",
    ) as output:
        output.write(f"{content}\n")


def main() -> None:
    plotting = True
    cli_args = _get_command_line_arguments()

    total_morfeus_sterimol_df = pd.DataFrame()

    for calculation_directory in pathlib.Path(
        cli_args.calculation_directory
    ).rglob("geom.xyz"):
        xyz_path = pathlib.Path(calculation_directory)
        r_group_name = pathlib.PurePath(xyz_path).parts[-3]
        substituent_name = pathlib.PurePath(xyz_path).parts[-2]
        print(xyz_path)
        calculator = smores.calculators.SterimolCalculator.init_from_file(
            molecule_path=xyz_path
        )
        calculator.set_core_smiles(
            core_smiles=pathlib.PurePath(xyz_path).parts[-3]
        )
        try:
            sterimol_parameters = calculator.calculate_sterimol_parameters(
                method_types=["morfeus"],
                radii_types=[
                    "alvarez",
                    "bondi",
                    "crc",
                    "rahm",
                    "pyykko",
                    "truhlar",
                    "streusel",
                ],
            )

            output_directory = pathlib.Path(cli_args.output_directory)
            write_sterimol_to_output(
                sterimol_parameters,
                output_directory,
                r_group_name,
                substituent_name,
            )
        except numpy.linalg.LinAlgError:
            print(f"{r_group_name}  {substituent_name}")

        if plotting:
            morfeus_sterimol_df = get_morfeus_sterimol_dataframe(
                sterimol_parameters, substituent_name, r_group_name
            )

            total_morfeus_sterimol_df = pd.concat(
                [total_morfeus_sterimol_df, morfeus_sterimol_df]
            )

    if plotting:
        total_morfeus_sterimol_df.columns = [
            "substituent",
            "r_group",
            "radii_type",
            "sterimol parameter",
            "value",
        ]
        plot = sns.stripplot(
            data=total_morfeus_sterimol_df,
            x="sterimol parameter",
            y="value",
            hue="radii_type",
            # dodge=True,
            # jitter=False,
        )
        fig = plot.get_figure()
        fig.savefig("morfeus_sterimol_xyz_performance.png")

        total_morfeus_sterimol_df.to_csv(
            "morfeus_xyz_sterimol_parameters.csv", index=False
        )


if __name__ == "__main__":
    main()
