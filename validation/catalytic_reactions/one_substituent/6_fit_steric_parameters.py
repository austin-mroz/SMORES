#!python
import argparse
import pathlib
import typing
from dataclasses import dataclass
import pandas as pd
import seaborn as sns


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

    smores_results = pd.read_csv(args.smores_results)
    experimental_results = pd.read_csv(args.experimental_results)
    for reaction in experimental_results["reaction"].unique():
        for sterimol_parameter_fit in _fit_sterimol_parameters(
            smores_results=smores_results,
            experimental_results=experimental_results,
            reaction=reaction,
        ):
            _plot_results(sterimol_parameter_fit)


def _plot_results(sterimol_parameter_fit: SterimolFit) -> None:
    sterimol_parameter_dataframe = pd.DataFrame({
        'experimental_ddGs': sterimol_parameter_fit.experimental_ddGs,
        'predicted_ddGs': sterimol_parameter_fit.predicted_ddGs,
    })
    sterimol_parameter_dataframe["substituent"] = sterimol_parameter_dataframe.index

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
        sterimol_parameter_dataframe["predicted_ddGs"].max(),
        
    )

    fig = plot.get_figure()
    fig.savefig(
        f"{radii_type}_fit_desymmetrization_results.png",
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
        experimental_results_reaction_subset = experimental_results[experimental_results["reaction"] == reaction]
        core = experimental_results_reaction_subset["core"].unique()[0]
        smores_results = smores_results[smores_results["core"] == core]

        experimental_ddG = experimental_results_reaction_subset[["ddG"]]
        sterimol_parameters = smores_results[parameter_combination]
        ols_fit = sm.OLS(experimental_ddG, sterimol_parameters).fit()
        print(ols_fit.summary())
        predicted_ddGs = _get_predicted_ddGs(sterimol_parameters, ols_fit)
        yield SterimolFit(
                name=reaction,
                core=core,
                L_coefficient = ols_fit.params.get("L"),
                B1_coefficient = ols_fit.params.get("B1"),
                B5_coefficient = ols_fit.params.get("B5"),
                fit_summary = ols_fit.summary(),
                r_squared = ols_fit.rsquared,
                experimental_ddGs = _get_experimental_ddGs(experimental_results_reaction_subset),
                predicted_ddGs = _get_predicted_ddGs(sterimol_parameters, ols_fit)
        )


def _get_experimental_ddGs(experimental_results_reaction_subset: pd.DataFrame) -> dict[str, float]:
    experimental_ddGs = {}
    for index, row in experimental_results_reaction_subset.iterrows():
        experimental_ddGs[row["substituent"].values[0]] = row["ddG"].values[0]
    return experimental_ddGs


def _get_predicted_ddGs(sterimol_parameters: pd.DataFrame, ols_fit: pd.core.series.Series) -> dict[str, float]:
    predicted_ddGs = {}
    for index, row in sterimol_parameters.iterrows():
        predicted_ddG = (
                (ols_fit.params.get("L", 0) * substituent["L"])
                + (ols_fit.params.get("B1", 0) * substituent["B1"])
                + (ols_fit.params.get("B5", 0) * substituent["B5"])
                )
        predicted_ddGs[row["substituent"].values[0]] = predicted_ddG
    return predicted_ddGs

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
        default=pathlib.Path.cwd() / "exprimental_ddG.csv"
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
