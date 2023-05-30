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

    df = _convert_database_to_dataframe(atomlite.Database(args.database))
    _plot_morfeus_sterimol_xyz_performance(df, args.output_directory)

    radii_types_to_exclude = ["crc", "truhlar", "alvarez", "rahm"]
    df = df[~df["radii_type"].isin(radii_types_to_exclude)]

    for core, core_df in df.groupby("core"):
        _plot_sterimol_parameter_by_radii_type(
            df=core_df,
            parameter="L",
            figure_path=args.output_directory / f"{core}_L.png",
        )
        _plot_sterimol_parameter_by_radii_type(
            df=core_df,
            parameter="B1",
            figure_path=args.output_directory / f"{core}_B1.png",
        )
        _plot_sterimol_parameter_by_radii_type(
            df=core_df,
            parameter="B5",
            figure_path=args.output_directory / f"{core}_B5.png",
        )


def _convert_database_to_dataframe(
    database: atomlite.Database,
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
                    radii_type,
                    entry.properties[f"{radii_type}_L"],
                    entry.properties[f"{radii_type}_B1"],
                    entry.properties[f"{radii_type}_B5"],
                ]
            ]
        )
        tempdf = pd.concat([df, nrow])
        df = pd.DataFrame(tempdf)

    df.columns = ["name", "core", "substituent", "radii_type", "L", "B1", "B5"]

    return df.reset_index()


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
    )
    sns.move_legend(
        plot,
        "upper left",
        bbox_to_anchor=(1.0, 1.0),
        title="substituent",
    )
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
        default=_get_output_directory() / "5_output",
    )
    return parser.parse_args()


def _get_output_directory() -> pathlib.Path:
    return pathlib.Path(str(pathlib.Path.cwd()).replace("work", "data"))


if __name__ == "__main__":
    main()
