import argparse
import logging
import pathlib
import smores
import streusel
import morfeus
import dbstep.Dbstep as db
import seaborn as sns
import pandas as pd
import glob


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("calculation_directory")
    return parser.parse_args()


def gen_DBStep_xyz_parameters(
        xyz_path: pathlib.Path,
) -> db.dbstep:

    dbstep_sterimol = db.dbstep(
            str(xyz_path),
            atom1=1,
            atom2=2,
            sterimol=True,
            measure='classic',
    )

    return dbstep_sterimol


def gen_morfeus_xyz_parameters(
        xyz_path: pathlib.Path,
        radii_type: str,
) -> morfeus.Sterimol:

    xyz_data = smores.utilities.read_xyz(xyz_path)

    streusel_radii = smores.utilities.get_streusel_radii(xyz_data)

    if radii_type == 'streusel':
        morfeus_sterimol = morfeus.Sterimol(
                xyz_data.elements,
                xyz_data.coordinates,
                1,
                2,
                radii=streusel_radii,
        )
    else:
        morfeus_sterimol = morfeus.Sterimol(
                xyz_data.elements,
                xyz_data.coordinates,
                1,
                2,
                radii_type=radii_type,
        )

    return morfeus_sterimol


def get_morfeus_sterimol_dataframe(
        xyz_path: pathlib.Path,
        substituent_name: str,
        r_group_name: str,
) -> pd.DataFrame:

    morfeus_radii_types = ['alvarez', 'bondi', 'crc', 'rahm', 'pyykko', 'truhlar', 'streusel']

    morfeus_sterimol_df = pd.DataFrame()

    for radii_type in morfeus_radii_types:
        morfeus_sterimol = gen_morfeus_xyz_parameters(xyz_path, radii_type)

        nrowLval = pd.DataFrame([[substituent_name, r_group_name, radii_type, 'L_value', morfeus_sterimol.L_value]])
        nrowB1val = pd.DataFrame([[substituent_name, r_group_name, radii_type, 'B1_value', morfeus_sterimol.B_1_value]])
        nrowB5val = pd.DataFrame([[substituent_name, r_group_name, radii_type, 'B5_value', morfeus_sterimol.B_5_value]])
        tempdf = pd.concat([morfeus_sterimol_df, nrowLval, nrowB1val, nrowB5val])
        morfeus_sterimol_df = pd.DataFrame(tempdf)
    morfeus_sterimol_df.columns = ['substituent', 'r_group', 'radii_type', 'sterimol parameter', 'value']

    return morfeus_sterimol_df


def main() -> None:
    cli_args = _get_command_line_arguments()

    total_morfeus_sterimol_df = pd.DataFrame()

    for calculation_directory in pathlib.Path(cli_args.calculation_directory).rglob('geom.xyz'):
        xyz_path = pathlib.Path(calculation_directory)
        r_group_name = pathlib.PurePath(xyz_path).parts[-3]
        substituent_name = pathlib.PurePath(xyz_path).parts[-2]

        morfeus_sterimol_df = get_morfeus_sterimol_dataframe(xyz_path, substituent_name, r_group_name)

        total_morfeus_sterimol_df = pd.concat([
            total_morfeus_sterimol_df,
            morfeus_sterimol_df
            ])

    total_morfeus_sterimol_df.columns = ['substituent', 'r_group', 'radii_type', 'sterimol parameter', 'value']
    plot = sns.stripplot(
            data=total_morfeus_sterimol_df,
            x="sterimol parameter",
            y="value",
            hue="radii_type",
            #dodge=True,
            #jitter=False,
    )
    fig = plot.get_figure()
    fig.savefig('morfeus_sterimol_xyz_performance.png')

    total_morfeus_sterimol_df.to_csv('morfeus_xyz_sterimol_parameters.csv', index=False)


if __name__ == '__main__':
    main()
