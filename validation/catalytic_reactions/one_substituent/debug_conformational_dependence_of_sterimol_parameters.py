#!python
import argparse
import json
import pathlib
import typing

import morfeus
import pandas as pd
import rdkit.Chem.AllChem as rdkit
import seaborn as sns


def main() -> None:
    args = _get_command_line_arguments()
    args.output_directory.mkdir(exist_ok=True, parents=True)
    conformer_directory = args.output_directory / "conformers"
    conformer_directory.mkdir(exist_ok=True, parents=True)

    excluded_atoms = None
    if args.substructure_indices is not None:
        with open(args.substructure_indices) as f:
            excluded_atoms = [
                index + 1 for index in json.load(f)["core_indices"]
            ]

    molecule = rdkit.MolFromMolFile(
        str(args.input_file),
        removeHs=False,
    )
    etkdg = rdkit.ETKDGv3()
    etkdg.randomSeed = args.random_seed
    rdkit.EmbedMultipleConfs(molecule, args.num_conformers, etkdg)

    Ls = []
    B1s = []
    B5s = []
    conformer_ids = []

    for conformer in molecule.GetConformers():
        sterimol = morfeus.Sterimol(
            elements=[atom.GetAtomicNum() for atom in molecule.GetAtoms()],
            coordinates=conformer.GetPositions(),
            dummy_index=args.dummy_index + 1,
            attached_index=args.attached_index + 1,
            excluded_atoms=excluded_atoms,
        )
        Ls.append(sterimol.L_value)
        B1s.append(sterimol.B_1_value)
        B5s.append(sterimol.B_5_value)
        conformer_id = conformer.GetId()
        conformer_ids.append(conformer_id)
        rdkit.MolToMolFile(
            molecule,
            str(
                conformer_directory / f"{conformer_id}.mol",
            ),
            confId=conformer_id,
            forceV3000=True,
        )
        rdkit.MolToXYZFile(
            molecule,
            str(
                conformer_directory / f"{conformer_id}.xyz",
            ),
            confId=conformer_id,
        )

    data_frame = pd.DataFrame(
        {
            "conformer_id": conformer_ids,
            "L": Ls,
            "B1": B1s,
            "B5": B5s,
        }
    )
    data_frame.to_csv(
        args.output_directory / "steric_parameters.csv",
        index=False,
    )
    plot_histogram(
        data_frame=data_frame,
        parameter="L",
        output_directory=args.output_directory,
    )
    plot_histogram(
        data_frame=data_frame,
        parameter="B1",
        output_directory=args.output_directory,
    )
    plot_histogram(
        data_frame=data_frame,
        parameter="B5",
        output_directory=args.output_directory,
    )


def plot_histogram(
    data_frame: pd.DataFrame,
    parameter: typing.Literal["L", "B1", "B5"],
    output_directory: pathlib.Path,
) -> None:

    plot = sns.histplot(
        data=data_frame,
        x=parameter,
    )
    plot.get_figure().savefig(
        output_directory / f"{parameter}.png",
        bbox_inches="tight",
    )
    plot.get_figure().clf()


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot sterimol parameters across ETKDG conformations."
    )
    parser.add_argument(
        "input_file",
        help="A MDL MOL file of a molecule.",
        type=pathlib.Path,
    )
    parser.add_argument(
        "dummy_index",
        help="The index of the dummy atom.",
        type=int,
    )
    parser.add_argument(
        "attached_index",
        help="The index of the attached atoms of the substituent.",
        type=int,
    )
    parser.add_argument(
        "output_directory",
        help="The directory in which the plots are saved.",
        type=pathlib.Path,
    )
    parser.add_argument(
        "--random_seed",
        help="The random seed to use for ETKDG.",
        type=int,
        default=4,
    )
    parser.add_argument(
        "--num_conformers",
        help="The number of conformers to generate.",
        type=int,
        default=500,
    )
    parser.add_argument(
        "--substructure_indices",
        help=(
            "A JSON file holding the indices of core and substituent atoms. "
            "If provided, the core atoms will be excluded from the sterimol "
            "parameter calculation."
        ),
        type=pathlib.Path,
        default=None,
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
