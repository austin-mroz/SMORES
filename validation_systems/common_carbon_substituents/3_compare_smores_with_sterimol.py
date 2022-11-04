import argparse
import enum
import pathlib
import textwrap

import morfeus


def main() -> None:
    args = _get_command_line_arguments()

    num_successes = 0
    num_failures = 0
    for xyz_file in args.xyz_files:
        sterimol = _get_sterimol(xyz_file)
        smores = _get_smores(xyz_file)
        result = _get_result(smores, sterimol)
        output = textwrap.dedent(
            f"""
            FILE: {xyz_file}
                SMORES L: {smores.L_value} STERIMOL L: {sterimol.L_value}
                SMORES B1: {smores.B_1_value} STERIMOL B1: {sterimol.B_1_value}
                SMORES B5: {smores.B_5_value} STERIMOL B5: {sterimol.B_5_value}
                RESULT: {result.value})
            """
        )
        print(output)
        if result == Result.SUCCESS:
            num_successes += 1
        else:
            num_failures += 1


class Result(enum.Enum):
    SUCCESS = "SUCCESS"
    FAILURE = "FAILURE"


def _get_sterimol(xyz_file: pathlib.Path) -> morfeus.Sterimol:
    elements, coordinates = morfeus.read_xyz(xyz_file)
    return morfeus.Sterimol(element, coordinates)


def _get_smores(xyz_file: pathlib.Path) -> morfeus.Sterimol:
    elements, coordinates = morfeus.read_xyz(xyz_file)
    return morfeus.Sterimol(element, coordinates)


def _get_result(
    smores: morfeus.Sterimol, sterimol: morfeus.Sterimol
) -> Result:
    pass


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Calculate and compare SMORES and sterimol parameters.",
    )

    default_input_directory = pathlib.Path(__file__).parent / "2_output"
    parser.add_argument(
        "--xyz_files",
        help="",
        nargs="+",
        type=pathlib.Path,
        default=default_input_directory.glob("**/*.xyz"),
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
