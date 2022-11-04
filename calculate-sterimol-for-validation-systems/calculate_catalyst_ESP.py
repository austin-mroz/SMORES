import argparse
import glob
import logging
import pathlib
import smores
import streusel


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("calculation_directory")
    return parser.parse_args()


def main() -> None:
    cli_args = _get_command_line_arguments()

    logging.basicConfig(
            filename='calculate_catalyst_ESP.log',
            level=logging.DEBUG,
    )

    for catalyst in pathlib.Path(cli_args.calculation_directory).rglob('xtbopt.xyz'):
        print(catalyst.resolve().parent)
        mol = smores.Molecule.init_from_xyz(catalyst)
        mol.calculate_electrostatic_potential(
                catalyst.resolve().parent,
        )


if __name__ == '__main__':
    main()
