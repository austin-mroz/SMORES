import argparse
import itertools
import pathlib

import pandas as pd
import seaborn as sns


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv")
    return parser.parse_args()


def plot(df: pd.DataFrame, output_directory: pathlib.Path) -> None:
    print(df)
    plot = sns.stripplot(
        data=df,
        x="radii_type",
        y="value",
        hue="substituent",
    )
    sns.move_legend(
        plot, "best", bbox_to_anchor=(1.25, 1.0), title="substituent"
    )
    fig = plot.get_figure()
    # fig.legend(loc="upper left")
    fig.savefig(str(output_directory), bbox_inches="tight")
    fig.clf()


def main() -> None:
    cli_args = _get_command_line_arguments()
    output_directory = pathlib.Path(cli_args.csv).parents[0]

    tot_df = pd.read_csv(cli_args.csv)

    l = tot_df[tot_df["sterimol parameter"] == "L_value"]
    b1 = tot_df[tot_df["sterimol parameter"] == "B1_value"]
    b5 = tot_df[tot_df["sterimol parameter"] == "B5_value"]

    for sterimol_params in zip([l, b1, b5], ["Bn_L", "Bn_B1", "Bn_B5"]):
        plot(
            sterimol_params[0],
            output_directory.joinpath(f"{sterimol_params[1]}.png"),
        )


if __name__ == "__main__":
    main()
