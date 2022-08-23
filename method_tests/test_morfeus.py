import argparse

import morfeus


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("xyz_file")
    return parser.parse_args()


def main() -> None:
    cli_args = _get_command_line_arguments()
    elements, coordinates = morfeus.read_xyz(cli_args.xyz_file)
    sterimol = morfeus.Sterimol(elements, coordinates, 2, 5)
    print(
        f"L: {sterimol.L} "
        f"-- B1: {sterimol.B_1_value} "
        f"-- B5: {sterimol.B_5_value}"
    )


if __name__ == "__main__":
    main()
