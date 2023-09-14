#!python

import argparse
import itertools
import pathlib

import atomlite
import pandas as pd
import seaborn as sns


def main() -> None:
    args = _get_command_line_arguments()
    args.output_directory.mkdir(parents=True, exist_ok=True)

    dipole_dict = _extract_substituent_dipole_moments(
        atomlite.Database(args.database)
    )

    df = _convert_database_to_dataframe(
        atomlite.Database(args.database),
        dipole_dict,
    )
    _plot_morfeus_sterimol_xyz_performance(df, args.output_directory)

    radii_types_to_exclude = ["crc", "truhlar", "alvarez", "rahm"]
    df = df[~df["radii_type"].isin(radii_types_to_exclude)]
    output_directory = args.output_directory / "main_text_fig2_radii"
    output_directory.mkdir(parents=True, exist_ok=True)
    for core, core_df in df.groupby("core"):
        _plot_sterimol_parameter_by_radii_type(
            df=core_df,
            parameter="L",
            figure_path=output_directory / f"{core}_L.png",
        )
        _plot_sterimol_parameter_by_radii_type(
            df=core_df,
            parameter="B1",
            figure_path=output_directory / f"{core}_B1.png",
        )
        _plot_sterimol_parameter_by_radii_type(
            df=core_df,
            parameter="B5",
            figure_path=output_directory / f"{core}_B5.png",
        )
        _plot_sterimol_parameter_by_dipole(
            df=core_df,
            parameter="L",
            figure_path=output_directory / f"{core}_dipole_c_by_subst_L.png",
        )
        _plot_sterimol_parameter_by_dipole(
            df=core_df,
            parameter="B1",
            figure_path=output_directory / f"{core}_dipole_c_by_subst_B1.png",
        )
        _plot_sterimol_parameter_by_dipole(
            df=core_df,
            parameter="B5",
            figure_path=output_directory / f"{core}_dipole_c_by_subst_B5.png",
        )
        _plot_sterimol_parameter_by_difference(
            df=core_df,
            parameter="L",
            figure_path=output_directory / f"{core}_difference.png",
        )


def _get_bondi_diff(
    df: pd.DataFrame,
    parameter: str,
) -> pd.DataFrame:
    bondi_df = df.loc[df["radii_type"] == "bondi"].reset_index()
    df = df[~df["radii_type"].isin(["bondi"])].reset_index()

    bondi_diff = []

    for index, row in df.iterrows():
        bondi_param = bondi_df.loc[
            bondi_df["substituent"] == row["substituent"], parameter
        ].iloc[0]
        bondi_diff.append(
            (abs(row[parameter] - bondi_param) / bondi_param) * 100
        )
    ncol = pd.DataFrame({"bondi_diff": bondi_diff})
    return pd.concat((df, ncol), axis=1)


def _plot_sterimol_parameter_by_difference(
    df: pd.DataFrame,
    parameter: str,
    figure_path: pathlib.Path,
) -> None:
    tdf_L = _get_bondi_diff(df, "L")
    tdf_B1 = _get_bondi_diff(df, "B1")
    tdf_B5 = _get_bondi_diff(df, "B5")

    plot = sns.stripplot(
        data=tdf_L,
        x="dipole",
        y="bondi_diff",
        hue="radii_type",
        edgecolor="black",
        linewidth=2,
        size=10,
        dodge=True,
        jitter=False,
    )

    sns.stripplot(
        data=tdf_B1,
        x="dipole",
        y="bondi_diff",
        hue="radii_type",
        marker="D",
        edgecolor="black",
        linewidth=2,
        size=10,
        dodge=True,
        jitter=False,
    )
    sns.stripplot(
        data=tdf_B5,
        x="dipole",
        y="bondi_diff",
        hue="radii_type",
        marker="H",
        edgecolor="black",
        linewidth=2,
        size=10,
        dodge=True,
        jitter=False,
    )

    sns.move_legend(
        plot,
        "upper left",
        bbox_to_anchor=(1.0, 1.0),
        title="substituent",
    )
    plot.get_figure().savefig(figure_path, bbox_inches="tight")
    plot.get_figure().clf()


def _extract_substituent_dipole_moments(
    database: atomlite.Database,
) -> dict[str, float]:
    dipole_dict = {}

    for entry in database.get_entries():
        print(entry.key)
        dipole_dict[entry.key] = _extract_dipole_from_xtb_output(
            pathlib.Path(entry.properties["esp_file"].removesuffix("ESP.cube"))
            / "output.dat"
        )
    return dipole_dict


def _extract_dipole_from_xtb_output(
    xtb_path: pathlib.Path,
) -> float:
    with open(xtb_path, "r") as xtb_output:
        lines = xtb_output.readlines()
    c = 0
    for line in lines:
        c += 1
        if "Electrostatic potential computed" in line:
            dipole_idx = c - 5
    return float(lines[dipole_idx].split()[2])


def _convert_database_to_dataframe(
    database: atomlite.Database,
    dipole_dict: dict[str, float],
) -> pd.DataFrame:
    df = pd.DataFrame()

    # get list of unique radii
    radii_types = [
        k.replace("_B1", "")
        for k, v in next(database.get_entries()).properties.items()
        if "_B1" in k
    ]

    for radii_type, entry in itertools.product(
        radii_types, database.get_entries()
    ):
        nrow = pd.DataFrame(
            [
                [
                    entry.key,
                    entry.properties["core"],
                    entry.properties["substituent"],
                    dipole_dict[entry.key],
                    radii_type,
                    entry.properties[f"{radii_type}_L"],
                    entry.properties[f"{radii_type}_B1"],
                    entry.properties[f"{radii_type}_B5"],
                ]
            ]
        )
        tempdf = pd.concat([df, nrow])
        df = pd.DataFrame(tempdf)

    df.columns = [
        "name",
        "core",
        "substituent",
        "dipole",
        "radii_type",
        "L",
        "B1",
        "B5",
    ]

    return df.reset_index()


def _plot_sterimol_parameter_by_dipole(
    df: pd.DataFrame,
    parameter: str,
    figure_path: pathlib.Path,
) -> None:
    plot = sns.stripplot(
        data=df,
        x="dipole",
        y=parameter,
        hue="substituent",
        edgecolor="black",
        linewidth=2,
        size=10,
        dodge=True,
        jitter=False,
    )
    sns.move_legend(
        plot,
        "upper left",
        bbox_to_anchor=(1.0, 1.0),
        title="substituent",
    )
    if parameter == "L":
        plot.set_ylim([0, 17])
    elif parameter == "B1":
        plot.set_ylim([0, 9])
    elif parameter == "B5":
        plot.set_ylim([0, 14])
    plot.get_figure().savefig(figure_path, bbox_inches="tight")
    plot.get_figure().clf()


def _plot_sterimol_parameter_by_radii_type(
    df: pd.DataFrame,
    parameter: str,
    figure_path: pathlib.Path,
) -> None:
    plot = sns.stripplot(
        data=df,
        x="radii_type",
        y=parameter,
        hue="substituent",
        edgecolor="black",
        linewidth=2,
        size=10,
        dodge=True,
        jitter=False,
    )
    sns.move_legend(
        plot,
        "upper left",
        bbox_to_anchor=(1.0, 1.0),
        title="substituent",
    )
    if parameter == "L":
        plot.set_ylim([0, 17])
    elif parameter == "B1":
        plot.set_ylim([0, 9])
    elif parameter == "B5":
        plot.set_ylim([0, 14])
    plot.get_figure().savefig(figure_path, bbox_inches="tight")
    plot.get_figure().clf()


def _plot_morfeus_sterimol_xyz_performance(
    df: pd.DataFrame,
    output_directory: pathlib.Path,
) -> None:
    plot = sns.stripplot(
        data=df.melt(
            id_vars=["radii_type"],
            value_vars=["L", "B1", "B5"],
            var_name="sterimol_parameter",
        ),
        x="sterimol_parameter",
        y="value",
        hue="radii_type",
        edgecolor="black",
        linewidth=2,
        size=10,
        dodge=True,
        jitter=True,
    )
    sns.move_legend(
        plot,
        "upper left",
        bbox_to_anchor=(1.0, 1.0),
        title="radii_type",
    )
    plot.get_figure().savefig(
        output_directory / "morfeus_sterimol_xyz_performance.png",
        bbox_inches="tight",
    )
    plot.get_figure().clf()


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d",
        "--database",
        help=(
            'An atomlite database file with properties: "core", "substituent", '
            '"xyz_file", "esp_file", dummy_index" and "attached_index".'
        ),
        type=pathlib.Path,
        default=pathlib.Path.cwd() / "common_carbon_substituents.db",
    )
    parser.add_argument(
        "-o",
        "--output_directory",
        help="The directory into which the results are written.",
        type=pathlib.Path,
        default=_get_output_directory() / "55_dipole_output",
    )
    return parser.parse_args()


def _get_output_directory() -> pathlib.Path:
    return pathlib.Path(str(pathlib.Path.cwd()).replace("work", "data"))


if __name__ == "__main__":
    main()
