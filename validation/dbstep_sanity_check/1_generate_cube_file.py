#!python
import argparse
import pathlib

import smores
import smores.psi4


def main() -> None:
    args = _get_command_line_arguments()
    ethane = smores.rdkit_from_smiles("CC")
    smores.psi4.calculate_electrostatic_potential(
        molecule=ethane,
        output_directory=args.output_directory,
        grid_origin=(-5, -5, -5),
        grid_length=10.0,
        num_voxels_per_dimension=50,
        optimize=True,
    )


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Creates a electrostatic potential cube file for ethane."
    )
    parser.add_argument(
        "-o",
        "--output_directory",
        type=pathlib.Path,
        default=pathlib.Path.cwd() / "1_output",
        help="The directory into which the output files are written.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
