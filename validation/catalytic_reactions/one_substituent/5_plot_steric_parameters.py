#!python

import argparse
import pathlib

import pandas as pd
import seaborn as sns


def main() -> None:
    args = _get_command_line_arguments()
    args.output_directory.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input_file)
    _plot_morfeus_sterimol_xyz_performance(df, args.output_directory)

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
    plot.get_figure().autofmt_xdate()
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
        "-i",
        "--input_file",
        help=(
            'A csv file with columns: "name", "core", '
            '"subsittuent", "radii_type", "L", "B1", "B5"'
        ),
        type=pathlib.Path,
        default=pathlib.Path.cwd() / "4_output" / "steric_parameters.csv",
    )
    parser.add_argument(
        "-o",
        "--output_directory",
        type=pathlib.Path,
        default=pathlib.Path.cwd() / "5_output",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
