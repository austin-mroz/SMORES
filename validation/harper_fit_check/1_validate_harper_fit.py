#!python
import argparse
import pathlib
import typing
from dataclasses import dataclass

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

    def add_row(
        self,
        steric_parameter_fit: StericParameterFit,
    ) -> "ResultsTable":
        self.L_coefficients.append(steric_parameter_fit.L_coefficient)
        self.B1_coefficients.append(steric_parameter_fit.B1_coefficient)
        self.B5_coefficients.append(steric_parameter_fit.B5_coefficient)
        self.intercepts.append(steric_parameter_fit.intercept)
        self.r_squareds.append(steric_parameter_fit.r_squared)
        return self


def main() -> None:
    args = _get_command_line_arguments()
    args.output_directory.mkdir(parents=True, exist_ok=True)

    literature_results = pd.read_csv(args.input_file)

    results_table = ResultsTable()

    for steric_parameter_fit in _fit_steric_parameters(literature_results):
        results_table.add_row(
            steric_parameter_fit=steric_parameter_fit,
        )
    results = pd.DataFrame(
        {
            "r_squared": results_table.r_squareds,
            "L_coefficient": results_table.L_coefficients,
            "B1_coefficient": results_table.B1_coefficients,
            "B5_coefficient": results_table.B5_coefficients,
            "intercept": results_table.intercepts,
        },
    ).sort_values(
        by=["r_squared"],
        ascending=False,
    )
    results.to_csv(
        args.output_directory / "r_squared.csv",
        index=False,
    )


def _powerset(*items: str) -> list[list[str]]:
    subsets: list[list[str]] = [[]]
    for item in items:
        subsets += [subset + [item] for subset in subsets]
    # Do not return the empty set.
    return subsets[1:]


def _fit_steric_parameters(
    data_frame: pd.DataFrame,
) -> typing.Iterator[StericParameterFit]:

    parameter_combinations = _powerset("L", "B1", "B5")
    for parameter_combination in parameter_combinations:
        X_unscaled = data_frame[parameter_combination]
        # X = (
        #     preprocessing.StandardScaler()
        #     .fit(X_unscaled)
        #     .transform(X_unscaled)
        # )
        X = X_unscaled
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


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fit sterimol parameters to experimental results.",
    )
    parser.add_argument(
        "--input_file",
        help="A csv file containing experimental ddG values.",
        type=pathlib.Path,
        default=pathlib.Path.cwd() / "desymm_nat_chem_sterimol_parameters.csv",
    )
    parser.add_argument(
        "-o",
        "--output_directory",
        type=pathlib.Path,
        default=pathlib.Path.cwd() / "1_output",
    )

    return parser.parse_args()


if __name__ == "__main__":
    main()
