import argparse
import glob
import logging
import pathlib
import itertools

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
    parser.add_argument("atom_idx_csv")
    parser.add_argument("structure_directory")
    parser.add_argument("output_directory")
    return parser.parse_args()


def get_morfeus_sterimol_dataframe(
    sterimol_parameters: list[SterimolParameter],
    substituent_name: str,
    r_group_name: str,
) -> pd.DataFrame:

    morfeus_sterimol_df = pd.DataFrame()

    for sterimol_ in list(itertools.chain(*sterimol_parameters)):

        nrowLval = pd.DataFrame(
            [
                [
                    sterimol_.method,
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
                    sterimol_.method,
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
                    sterimol_.method,
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
        "method",
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
    for sterimol_ in list(itertools.chain(*sterimol_parameters)):
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


def get_atom_idx(atom_idx_df: pd.DataFrame, substituent: str,) -> dict[int, any]:
    substituent_atom_idx = atom_idx_df[atom_idx_df['substituent'] == substituent]
    return {'atom_1': substituent_atom_idx['atom_1'].values,
            'atom_2': substituent_atom_idx['atom_2'].values
            }


def main() -> None:
    plotting = True
    cli_args = _get_command_line_arguments()

    total_morfeus_sterimol_df = pd.DataFrame()

    atom_idx_df = pd.read_csv(cli_args.atom_idx_csv)

    for calculation_directory in pathlib.Path(
        cli_args.structure_directory
    ).rglob("geom.xyz"):
        xyz_path = pathlib.Path(calculation_directory)
        r_group_name = pathlib.PurePath(xyz_path).parts[-3]
        substituent_name = pathlib.PurePath(xyz_path).parts[-2]

        print(substituent_name)

        calculator = smores.calculators.SterimolCalculator.init_from_file(
            molecule_path=xyz_path
        )

        atom_idx_dict = get_atom_idx(atom_idx_df, substituent_name) 
        calculator.set_atom_idx(
            atom_1_idx=atom_idx_dict['atom_1'][0],
            atom_2_idx=atom_idx_dict['atom_2'][0],
        )
        try:
            sterimol_parameters = calculator.calculate_sterimol_parameters(
                method_types=["morfeus", "db"],
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
            "method",
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
        fig.savefig("sterimol_xyz_performance.png")

        total_morfeus_sterimol_df.to_csv(
            "sterimol_parameters.csv", index=False
        )


if __name__ == "__main__":
    main()
