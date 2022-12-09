#!python
import argparse
import csv
import json
import pathlib

import rdkit.Chem as rdkit

import smores


def main() -> None:
    args = _get_command_line_arguments()
    args.output_directory.mkdir(parents=True, exist_ok=True)

    cores_allylation = {
        smores.rdkit_from_smiles(
            "CC(C)[C@H](NC(=O)C1CCCN1C(=O)OBr)C3=N[C@H](Cc2ccccc2)CO3"
        ): "allylation_catalyst",
    }
    substituents_allylation = {
        smores.rdkit_from_smiles("BrC"): "Me",
        smores.rdkit_from_smiles("BrCC"): "Et",
        smores.rdkit_from_smiles("CC(C)Br"): "i-Pr",
        smores.rdkit_from_smiles("CC(C)(C)Br"): "tBu",
        smores.rdkit_from_smiles("CCCC(Br)CCC"): "CH(Pr)2",
        smores.rdkit_from_smiles("CCC(Br)(CC)CC"): "CEt3",
        smores.rdkit_from_smiles("CC(C)C(Br)C(C)C"): "CH(i-Pr)2",
        smores.rdkit_from_smiles("BrC23CC1CC(CC(C1)C2)C3"): "1-Ad",
    }

    cores_desymmetrization = {
        smores.rdkit_from_smiles(
            "Oc2ccc(C(Br)c1ccc(O)cc1)cc2"
        ): "desymmetrization_reactant",
    }
    substituents_desymmetrization = {
        smores.rdkit_from_smiles("BrC"): "Me",
        smores.rdkit_from_smiles("BrCC"): "Et",
        smores.rdkit_from_smiles("Brc1ccccc1"): "Ph",
        smores.rdkit_from_smiles("BrCc1ccccc1"): "Bn",
        smores.rdkit_from_smiles("N#CBr"): "Cy",
        smores.rdkit_from_smiles("BrC(C)(C)CC"): "CH2-tBu",
        smores.rdkit_from_smiles("CCC(Br)CC"): "CHEt2",
        smores.rdkit_from_smiles("BrC(C)CC"): "CH2-iPr",
        smores.rdkit_from_smiles("BrC(c1ccccc1)c2ccccc2"): "CH(Ph)2",
        smores.rdkit_from_smiles("CC(C)Br"): "i-Pr",
        smores.rdkit_from_smiles("CC(C)(C)Br"): "tBu",
        smores.rdkit_from_smiles("BrC23CC1CC(CC(C1)C2)C3"): "1-Ad",
    }

    _write_structures(
        reaction_name="allylation",
        output_directory=args.output_directory,
        cores=cores_allylation,
        substituents=substituents_allylation,
    )
    _write_structures(
        reaction_name="desymmetrization",
        output_directory=args.output_directory,
        cores=cores_desymmetrization,
        substituents=substituents_desymmetrization,
    )


def _write_structures(
    reaction_name: str,
    output_directory: pathlib.Path,
    cores: dict[rdkit.Mol, str],
    substituents: dict[rdkit.Mol, str],
) -> None:
    structures_directory = output_directory / reaction_name
    structures_directory.mkdir(parents=True, exist_ok=True)

    with open(structures_directory / "xyz_files.csv", "w") as csv_file:
        writer = csv.DictWriter(
            csv_file,
            fieldnames=[
                "reaction_name",
                "name",
                "core",
                "substituent",
                "smiles",
                "dummy_index",
                "attached_index",
                "xyz_file",
            ],
        )
        writer.writeheader()
        for combo in smores.combine(cores, substituents):
            reaction_name = reaction_name
            name = f"{cores[combo.core]}_{substituents[combo.substituent]}"
            xyz_file = structures_directory / f"{name}.xyz"
            rdkit.MolToXYZFile(combo.product, str(xyz_file))
            writer.writerow(
                {
                    "reaction_name": reaction_name,
                    "name": name,
                    "core": cores[combo.core],
                    "substituent": substituents[combo.substituent],
                    "smiles": rdkit.MolToSmiles(
                        rdkit.RemoveHs(combo.product),
                        canonical=True,
                    ),
                    "dummy_index": combo.dummy_index,
                    "attached_index": combo.attached_index,
                    "xyz_file": xyz_file.resolve(),
                },
            )
            with open(structures_directory / f"{name}.json", "w") as f:
                json.dump(
                    {
                        "core_indices": combo.core_indices,
                        "substituent_indices": combo.substituent_indices,
                    },
                    f,
                )


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate the structures of systems which are to be validated."
        ),
    )

    parser.add_argument(
        "-o",
        "--output_directory",
        default=pathlib.Path.cwd() / "1_output",
        type=pathlib.Path,
        help="The directory into which the output files are written.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
