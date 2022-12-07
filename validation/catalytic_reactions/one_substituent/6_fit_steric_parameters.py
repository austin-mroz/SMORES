#!python
import argparse
import pathlib
import typing
from dataclasses import dataclass

import pandas as pd
import seaborn as sns
from sklearn.linear_model import LinearRegression


@dataclass(frozen=True, slots=True)
class StericParameterFit:
    L_coefficient: float | None
    B1_coefficient: float | None
    B5_coefficient: float | None
    r_squared: float
    results: pd.DataFrame


class ResultsTable:
    def __init__(self) -> None:
        self.L_coefficients: list[float | None] = []
        self.B1_coefficients: list[float | None] = []
        self.B5_coefficients: list[float | None] = []
        self.r_squareds: list[float] = []
        self.reactions: list[str] = []
        self.radii_types: list[str] = []

    def add_row(
        self,
        steric_parameter_fit: StericParameterFit,
        reaction: str,
        radii_type: str,
    ) -> "ResultsTable":
        self.L_coefficients.append(steric_parameter_fit.L_coefficient)
        self.B1_coefficients.append(steric_parameter_fit.B1_coefficient)
        self.B5_coefficients.append(steric_parameter_fit.B5_coefficient)
        self.r_squareds.append(steric_parameter_fit.r_squared)
        self.reactions.append(reaction)
        self.radii_types.append(radii_type)
        return self


def main() -> None:
    args = _get_command_line_arguments()
    args.output_directory.mkdir(parents=True, exist_ok=True)

    steric_parameters = pd.read_csv(args.steric_parameters)
    experimental_ddGs = pd.read_csv(args.experimental_ddGs)

    results_table = ResultsTable()
    for radii_type in steric_parameters["radii_type"].unique():
        args.output_directory.joinpath(radii_type).mkdir(
            exist_ok=True,
            parents=True,
        )

        for reaction in experimental_ddGs["reaction"].unique():
            experimental_ddGs_of_reaction = experimental_ddGs[
                experimental_ddGs["reaction"] == reaction
            ]
            (core,) = experimental_ddGs_of_reaction["core"].unique()
            for steric_parameter_fit in _fit_steric_parameters(
                data_frame=pd.merge(
                    left=steric_parameters[
                        (steric_parameters["core"] == core)
                        & (steric_parameters["radii_type"] == radii_type)
                    ],
                    right=experimental_ddGs_of_reaction,
                    on=["substituent", "core"],
                    how="inner",
                ),
            ):
                _plot_results(
                    steric_parameter_fit=steric_parameter_fit,
                    output_directory=args.output_directory / radii_type,
                )
                results_table.add_row(
                    steric_parameter_fit=steric_parameter_fit,
                    reaction=reaction,
                    radii_type=radii_type,
                )

    results = pd.DataFrame(
        {
            "r_squared": results_table.r_squareds,
            "reaction": results_table.reactions,
            "radii_type": results_table.radii_types,
            "L_coefficient": results_table.L_coefficients,
            "B1_coefficient": results_table.B1_coefficients,
            "B5_coefficient": results_table.B5_coefficients,
        },
    ).sort_values(
        by=["reaction", "radii_type", "r_squared"],
        ascending=False,
    )
    results.to_csv(
        args.output_directory / "r_squared.csv",
        index=False,
    )

    max_r_squared_rows = (
        results.groupby(["reaction", "radii_type"])["r_squared"].transform(max)
        == results["r_squared"]
    )
    results[max_r_squared_rows].sort_values(
        by=["reaction", "r_squared"], ascending=False
    ).to_csv(
        args.output_directory / "r_squared_summary.csv",
        index=False,
    )


def _fit_steric_parameters(
    data_frame: pd.DataFrame,
) -> typing.Iterator[StericParameterFit]:

    parameter_combinations = _powerset("L", "B1", "B5")
    for parameter_combination in parameter_combinations:
        X = data_frame[parameter_combination]
        y = data_frame["ddG"]
        ols_fit = LinearRegression().fit(X, y)
        coefficients = dict(zip(parameter_combination, ols_fit.coef_))
        predicted_ddGs = pd.DataFrame(
            {
                "substituent": data_frame["substituent"],
                "predicted_ddG": ols_fit.predict(X),
            },
        )
        yield StericParameterFit(
            L_coefficient=coefficients.get("L"),
            B1_coefficient=coefficients.get("B1"),
            B5_coefficient=coefficients.get("B5"),
            r_squared=ols_fit.score(X, y),
            results=pd.merge(
                left=data_frame,
                right=predicted_ddGs,
                on="substituent",
                how="inner",
            ),
        )


def _plot_results(
    steric_parameter_fit: StericParameterFit,
    output_directory: pathlib.Path,
) -> None:

    plot = sns.scatterplot(
        data=steric_parameter_fit.results,
        x="ddG",
        y="predicted_ddG",
        hue="substituent",
    )
    sns.move_legend(
        plot,
        "upper left",
        bbox_to_anchor=(1.0, 1.0),
        title="substituent",
    )

    fit_equation = (
        f"{steric_parameter_fit.L_coefficient} L "
        f"+ {steric_parameter_fit.B1_coefficient} B1 "
        f"+ {steric_parameter_fit.B5_coefficient} B5"
    )

    plot.text(
        1.05,
        0.0,
        fit_equation,
        ha="left",
        va="bottom",
        transform=plot.transAxes,
    )
    plot.text(
        1.05,
        -0.05,
        f"R2 = {steric_parameter_fit.r_squared}",
        ha="left",
        va="bottom",
        transform=plot.transAxes,
    )
    (reaction,) = steric_parameter_fit.results["reaction"].unique()
    (core,) = steric_parameter_fit.results["core"].unique()
    basename = "_".join(
        map(
            str,
            [
                reaction,
                core,
                steric_parameter_fit.L_coefficient,
                steric_parameter_fit.B1_coefficient,
                steric_parameter_fit.B5_coefficient,
            ],
        )
    )
    output = output_directory / basename
    steric_parameter_fit.results.to_csv(
        f"{output}.csv",
        columns=[
            "reaction",
            "core",
            "substituent",
            "ddG",
            "predicted_ddG",
        ],
        index=False,
    )
    fig = plot.get_figure()
    fig.savefig(
        f"{output}.png",
        bbox_inches="tight",
    )
    fig.clf()


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
        "--steric_parameters",
        help=(
            "A csv file containing steric parameters used to "
            "predict ddG values."
        ),
        type=pathlib.Path,
        default=pathlib.Path.cwd() / "4_output" / "steric_parameters.csv",
    )
    parser.add_argument(
        "--experimental_ddGs",
        help="A csv file containing experimental ddG values.",
        type=pathlib.Path,
        default=pathlib.Path.cwd() / "experimental_ddG.csv",
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
