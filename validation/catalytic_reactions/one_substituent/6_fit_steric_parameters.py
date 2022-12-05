#!python
import argparse
import pathlib
import typing
from dataclasses import dataclass
import statsmodels.api as sm
import pandas as pd
import seaborn as sns
from sklearn.linear_model import LinearRegression


_OUTPUT_CSV_COLUMNS = (
    "name",
    "core",
    "substituent",
    "smiles",
    "xyz_file",
    "dummy_index",
    "attached_index",
    "radii_type",
    "L",
    "B1",
    "B5",
)


@dataclass(frozen=True, slots=True)
class SterimolFit:
    name: str
    core: str
    L_coefficient: float | None
    B1_coefficient: float | None
    B5_coefficient: float | None
    fit_summary: list[str]
    r_squared: float
    experimental_ddGs: dict[str, float]
    predicted_ddGs: dict[str, float]


def main() -> None:
    args = _get_command_line_arguments()
    args.output_directory.mkdir(parents=True, exist_ok=True)

    all_smores_results = pd.read_csv(args.smores_results)
    smores_results = all_smores_results[all_smores_results["radii_type"] == "bondi"]
    experimental_results = pd.read_csv(args.experimental_results)
    for reaction in experimental_results["reaction"].unique():
        for sterimol_parameter_fit in _fit_sterimol_parameters(
            smores_results=smores_results,
            experimental_results=experimental_results,
            reaction=reaction,
        ):
            _plot_results(sterimol_parameter_fit, args.output_directory)


def _plot_results(
    sterimol_parameter_fit: SterimolFit, output_directory: pathlib.Path
) -> None:
    sterimol_parameter_dataframe = pd.DataFrame(
        {
            "experimental_ddGs": sterimol_parameter_fit.experimental_ddGs,
            "predicted_ddGs": sterimol_parameter_fit.predicted_ddGs,
        }
    )
    sterimol_parameter_dataframe[
        "substituent"
    ] = sterimol_parameter_dataframe.index
    plot = sns.scatterplot(
        data=sterimol_parameter_dataframe,
        x="experimental_ddGs",
        y="predicted_ddGs",
        hue="substituent",
    )
    sns.move_legend(
        plot,
        "upper left",
        bbox_to_anchor=(1.0, 1.0),
        title="substituent",
    )

    fit_equation = f"{sterimol_parameter_fit.L_coefficient}L + {sterimol_parameter_fit.B1_coefficient}B1 + {sterimol_parameter_fit.B5_coefficient}B5"

    plot.text(
        0,
        0,# sterimol_parameter_dataframe["predicted_ddGs"].max(),
        fit_equation,
    )
    plot.text(
        0,
        0, # sterimol_parameter_dataframe["predicted_ddGs"].max() - 0.5,
        sterimol_parameter_fit.r_squared,
    )
    basename = "_".join(
            map(str,
        [
            sterimol_parameter_fit.name,
            sterimol_parameter_fit.core,
            sterimol_parameter_fit.L_coefficient,
            sterimol_parameter_fit.B1_coefficient,
            sterimol_parameter_fit.B5_coefficient,
        ])
    )
    fig = plot.get_figure()
    fig.savefig(
        str(output_directory / f"{basename}.png"),
        bbox_inches="tight",
    )
    fig.clf()


def _fit_sterimol_parameters(
    smores_results: pd.DataFrame,
    experimental_results: pd.DataFrame,
    reaction: str,
) -> typing.Iterator[SterimolFit]:
    parameter_combinations = _powerset("L", "B1", "B5")
    for parameter_combination in parameter_combinations:
        experimental_results_reaction_subset = experimental_results[
            experimental_results["reaction"] == reaction
        ]
        core = experimental_results_reaction_subset["core"].unique()[0]
        smores_results_reaction_subset = smores_results[smores_results["core"] == core]

        experimental_ddG = experimental_results_reaction_subset[["ddG"]]
        sterimol_parameters = smores_results_reaction_subset[parameter_combination]
        ols_fit = LinearRegression().fit(sterimol_parameters, experimental_ddG)
        ols_fit_params = ols_fit.get_params()
        yield SterimolFit(
            name=reaction,
            core=core,
            L_coefficient=ols_fit_params.get("L"),
            B1_coefficient=ols_fit_params.get("B1"),
            B5_coefficient=ols_fit_params.get("B5"),
            fit_summary= "summary", # ols_fit.summary(),
            r_squared=ols_fit.score(sterimol_parameters, experimental_ddG),
            experimental_ddGs=_get_experimental_ddGs(
                experimental_results_reaction_subset
            ),
            predicted_ddGs=_get_predicted_ddGs(sterimol_parameters, smores_results_reaction_subset, ols_fit)
        )


def _get_predicted_ddGs(
        sterimol_parameters: pd.DataFrame,
        smores_results_reaction_subset: pd.DataFrame,
        ols_fit: LinearRegression,
) -> dict[str, float]:
    predicted_ddGs = {}
    ols_fit_predicted_values = pd.DataFrame(ols_fit.predict(sterimol_parameters))
    smores_substituents = smores_results_reaction_subset["substituent"]
    return dict(zip(smores_substituents, ols_fit_predicted_values))


def _get_experimental_ddGs(
    experimental_results_reaction_subset: pd.DataFrame,
) -> dict[str, float]:
    experimental_ddGs = {}
    for index, row in experimental_results_reaction_subset.iterrows():
        experimental_ddGs[row["substituent"]] = row["ddG"]
    return experimental_ddGs


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
        "--smores_results",
        help="A csv file containing sterimol parameters.",
        type=pathlib.Path,
        default=pathlib.Path.cwd() / "4_output" / "steric_parameters.csv",
    )
    parser.add_argument(
        "--experimental_results",
        help="A csv file containing experimental ddG.",
        type=pathlib.Path,
        default=pathlib.Path.cwd() / "exprimental_ddG.csv",
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
