import argparse
import pathlib
import morfeus


def main() -> None:
    args = _get_command_line_arguments()

    for cube_file in args.cube_file:
        pass


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Calculate and compare SMORES and sterimol parameters.",
    )
    parser.add_argument(
        "--cube_files",
        nargs="+",
        type=pathlib.Path,
        default=pathlib.Path(__file__).parent,
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
