#!python
import argparse
import itertools
import pathlib
import typing
from dataclasses import dataclass

import atomlite
import pandas as pd
import seaborn as sns
from sklearn import preprocessing
from sklearn.linear_model import LinearRegression


@dataclass(frozen=True, slots=True)
class StericParameterFit:
    L_coefficient: float | None
    B1_coefficient: float | None
    B5_coefficient: float | None
    intercept: float
    r_squared: float
    results: pd.DataFrame


class ResultsTable:
    def __init__(self) -> None:
        self.L_coefficients: list[float | None] = []
        self.B1_coefficients: list[float | None] = []
        self.B5_coefficients: list[float | None] = []
        self.intercepts: list[float] = []
        self.r_squareds: list[float] = []
        self.radii_types: list[str] = []

    def add_row(
        self,
        steric_parameter_fit: StericParameterFit,
        radii_type: str,
    ) -> "ResultsTable":
        self.L_coefficients.append(steric_parameter_fit.L_coefficient)
        self.B1_coefficients.append(steric_parameter_fit.B1_coefficient)
        self.B5_coefficients.append(steric_parameter_fit.B5_coefficient)
        self.intercepts.append(steric_parameter_fit.intercept)
        self.r_squareds.append(steric_parameter_fit.r_squared)
        self.radii_types.append(radii_type)
        return self


def main() -> None:
    args = _get_command_line_arguments()
    args.output_directory.mkdir(parents=True, exist_ok=True)

    experimental_ddGs_of_reaction = pd.read_csv(args.experimental_results)

    steric_parameters = _convert_database_to_dataframe(
        atomlite.Database(args.database)
    )

    results_table = ResultsTable()
    for radii_type in steric_parameters["radii_type"].unique():
        args.output_directory.joinpath(radii_type).mkdir(
            exist_ok=True,
            parents=True,
        )

        for steric_parameter_fit in _fit_steric_parameters(
            data_frame=pd.merge(
                left=steric_parameters[
                    (steric_parameters["radii_type"] == radii_type)
                ],
                right=experimental_ddGs_of_reaction,
                on=["substituent"],
                how="inner",
            ),
        ):
            _plot_results(
                steric_parameter_fit=steric_parameter_fit,
                output_directory=args.output_directory / radii_type,
            )
            results_table.add_row(
                steric_parameter_fit=steric_parameter_fit,
                radii_type=radii_type,
            )

    results = pd.DataFrame(
        {
            "r_squared": results_table.r_squareds,
            "radii_type": results_table.radii_types,
            "L_coefficient": results_table.L_coefficients,
            "B1_coefficient": results_table.B1_coefficients,
            "B5_coefficient": results_table.B5_coefficients,
            "intercept": results_table.intercepts,
        },
    ).sort_values(
        by=["radii_type", "r_squared"],
        ascending=False,
    )
    results.to_csv(
        args.output_directory / "r_squared.csv",
        index=False,
    )

    max_r_squared_rows = (
        results.groupby(["radii_type"])["r_squared"].transform(max)
        == results["r_squared"]
    )
    results[max_r_squared_rows].sort_values(
        by=["r_squared"], ascending=False
    ).to_csv(
        args.output_directory / "r_squared_summary.csv",
        index=False,
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


def _fit_steric_parameters(
    data_frame: pd.DataFrame,
) -> typing.Iterator[StericParameterFit]:
    parameter_combinations = _powerset("L", "B1", "B5")
    for parameter_combination in parameter_combinations:
        X_unscaled = data_frame[parameter_combination]
        X = (
            preprocessing.StandardScaler()
            .fit(X_unscaled)
            .transform(X_unscaled)
        )
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
            intercept=ols_fit.intercept_,
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
    basename = "_".join(
        map(
            str,
            [
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
        "-d",
        "--database",
        help=(
            'An atomlite database file with properties: "core", "substituent", '
            '"xyz_file", "esp_file", dummy_index" and "attached_index".'
        ),
        type=pathlib.Path,
        default=pathlib.Path.cwd() / "catalyst_systems.db",
    )
    parser.add_argument(
        "-e",
        "--experimental_results",
        help=('A csv file with columns: "key", "substituent", and' '"yield".'),
        type=pathlib.Path,
        default=pathlib.Path.cwd() / "experimental_ddGs.csv",
    )
    parser.add_argument(
        "-o",
        "--output_directory",
        help="The directory into which the results are written.",
        type=pathlib.Path,
        default=_get_output_directory() / "6_output",
    )
    return parser.parse_args()


def _get_output_directory() -> pathlib.Path:
    return pathlib.Path(str(pathlib.Path.cwd()).replace("work", "data"))


if __name__ == "__main__":
    main()
