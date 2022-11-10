import argparse
import numpy as np
import pandas as pd
import pathlib
import seaborn as sns
import statsmodels.api as sm


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("sterimol_parameters")
    parser.add_argument("experimental_results")
    parser.add_argument("output_name")
    return parser.parse_args()


def configure_sterimol_parameter_dataframe(
        sterimol_parameters: pd.DataFrame,
        experimental_results: pd.DataFrame,
        method: str,
        radii_type: str,
) -> pd.DataFrame:

    sterimol_parameters = sterimol_parameters[
            sterimol_parameters["method"] == method]
    sterimol_parameters = sterimol_parameters[
            sterimol_parameters["radii_type"] == radii_type]

    configured_sterimol_parameters = pd.DataFrame()

    for index, row in experimental_results.iterrows():
        substituent_sterimol_parameters = sterimol_parameters[
                sterimol_parameters["substituent"] == row["substituent"]]

        if not substituent_sterimol_parameters.empty:
            [L, B1, B5] = _get_sterimol_parameter_values(substituent_sterimol_parameters)

            nrow = pd.DataFrame([[
                substituent_sterimol_parameters["method"].values[0],
                substituent_sterimol_parameters["substituent"].values[0],
                substituent_sterimol_parameters["radii_type"].values[0],
                L,
                B1,
                B5,
                row["ddG"]
            ]])
            tempdf = pd.concat([configured_sterimol_parameters, nrow])
            configured_sterimol_parameters = pd.DataFrame(tempdf)
    configured_sterimol_parameters.columns = ['method', 'substituent', 'radii_type', 'L_value', 'B1_value', 'B5_value', 'ddG']

    return configured_sterimol_parameters


def _get_sterimol_parameter_values(
    sterimol_parameter_df: pd.DataFrame,
) -> list[float]:
    parameters = []
    for parameter in ["L_value", "B1_value", "B5_value"]:
        parameters.append(
                sterimol_parameter_df[
                    sterimol_parameter_df["sterimol parameter"] == parameter]["value"].values[0]
        )
    return parameters


def get_stats(X: pd.DataFrame, y: pd.DataFrame) -> None:
    results = sm.OLS(y, X).fit()
    print(results.summary())
    return results


def _get_predicted_ddG(
        calculated_parameters: pd.DataFrame,
        sterimol_coefficients: pd.core.series.Series,
) -> float:
    calculated_B1_value = calculated_parameters["B1_value"]

    B1_coeff = sterimol_coefficients["B1_value"]

    predicted_ddG = (B1_coeff * calculated_B1_value) 

    return predicted_ddG


def add_predicted_ddGs(
        sterimol_coefficients: pd.core.series.Series,
        sterimol_results: pd.DataFrame,
) -> pd.DataFrame:
    predicted_ddGs = []
    for index, row in sterimol_results.iterrows():
        predicted_ddG = _get_predicted_ddG(row, sterimol_coefficients)
        predicted_ddGs.append(predicted_ddG)
    
    sterimol_results['predicted_ddG'] = pd.Series(predicted_ddGs).values

    return sterimol_results


def plot_scatter_stats(
        sterimol_results: pd.DataFrame,
        radii_type: str,
        output_name: str,
) -> None:
    plot = sns.scatterplot(
            data=sterimol_results,
            x="ddG",
            y="predicted_ddG",
            hue="substituent",
    )
    plot.set(ylim=(0,2.0))
    plot.set(xlim=(0,1.8))
    fig = plot.get_figure()
    fig.savefig(f"{output_name}_{radii_type}_fit_desymmetrization_results.png")
    fig.clf()


def plot_regplot_stats(
        sterimol_results: pd.DataFrame,
        output_name:str,
) -> None:
    plot = sns.lmplot(
            data=sterimol_results,
            x="ddG",
            y="predicted_ddG",
            hue="radii_type",
    )
    plot.set(ylim=(0, 2.0))
    plot.set(xlim=(0.3, 1.8))
    # fig = plot.get_figure()
    plot.savefig(f"{output_name}_regplot_fit_desymmetrization_results.png")
    # plot.clf()


def main() -> None:
    cli_args = _get_command_line_arguments()

    sterimol_parameters = pd.read_csv(cli_args.sterimol_parameters)
    experimental_results = pd.read_csv(cli_args.experimental_results)

    radii_types = [
            "alvarez",
            "bondi",
            "crc",
            "rahm",
            "pyykko",
            "truhlar",
            "streusel"
    ]

    all_sterimol_results = pd.DataFrame()

    for radii_type in radii_types:
        sterimol_parameters_for_fitting = configure_sterimol_parameter_dataframe(
                sterimol_parameters,
                experimental_results,
                "morfeus",
                radii_type,
        )
        X = sterimol_parameters_for_fitting[[
            "B1_value",
            ]]
        y = sterimol_parameters_for_fitting[["ddG"]]

        sterimol_coefficients = get_stats(X, y)
        sterimol_results = add_predicted_ddGs(sterimol_coefficients.params, sterimol_parameters_for_fitting)

        plot_scatter_stats(sterimol_results, radii_type, cli_args.output_name)

        all_sterimol_results = all_sterimol_results.append(sterimol_results)
    plot_regplot_stats(all_sterimol_results, cli_args.output_name)


if __name__ == "__main__":
    main()
