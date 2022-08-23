import argparse

from streusel import gaussian_cube


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("cube_file")
    return parser.parse_args()


def main() -> None:
    cli_args = _get_command_line_arguments()
    mol = gaussian_cube.Molecule(cli_args.cube_file)
    mol.get_efield()
    mol.sample_efield()
    print(f"VOLMUME {mol.vol}")


if __name__ == "__main__":
    main()
