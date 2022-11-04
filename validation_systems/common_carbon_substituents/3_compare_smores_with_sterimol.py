import argparse
import enum
import pathlib
import textwrap

import morfeus

import smores


def main() -> None:
    args = _get_command_line_arguments()

    num_successes = 0
    num_failures = 0
    for xyz_file in args.xyz_files:
        sterimol = _get_sterimol(xyz_file)
        smores = _get_smores(xyz_file)
        result = _get_result(smores, sterimol, args.success_tolerance)
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
    print(f"TOTAL SUCCESSES: {num_successes}")
    print(f"TOTAL FAILURES: {num_failures}")


class Result(enum.Enum):
    SUCCESS = "SUCCESS"
    FAILURE = "FAILURE"


def _get_sterimol(xyz_file: pathlib.Path) -> morfeus.Sterimol:
    elements, coordinates = morfeus.read_xyz(xyz_file)
    return morfeus.Sterimol(
        elements=elements,
        coordinates=coordinates,
        dummy_index=0,
        attached_index=1,
    )


def _get_smores(xyz_file: pathlib.Path) -> morfeus.Sterimol:
    elements, coordinates = morfeus.read_xyz(xyz_file)
    return morfeus.Sterimol(
        elements=elements,
        coordinates=coordinates,
        dummy_index=0,
        attached_index=1,
        radii=[smores.streusel_radii[element] for element in elements],
    )


def _get_result(
    smores: morfeus.Sterimol,
    sterimol: morfeus.Sterimol,
    success_tolerance: float,
) -> Result:

    l_diff = abs(smores.L_value - sterimol.L_value)
    b1_diff = abs(smores.B_1_value - sterimol.B_1_value)
    b5_diff = abs(smores.B_5_value - sterimol.B_5_value)
    if (
        l_diff < success_tolerance
        and b1_diff < success_tolerance
        and b5_diff < success_tolerance
    ):
        return Result.SUCCESS
    else:
        return Result.FAILURE


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Calculate and compare SMORES and sterimol parameters.",
    )

    default_input_directory = pathlib.Path(__file__).parent / "2_output"
    parser.add_argument(
        "--xyz_files",
        help=(
            "The xyz files for which SMORES and sterimol parameters "
            "are compared."
        ),
        nargs="+",
        type=pathlib.Path,
        default=default_input_directory.glob("**/*.xyz"),
    )
    parser.add_argument(
        "--success_tolerance",
        help=(
            "The maximum allowed difference between SMORES and sterimol "
            "parameters for the validation to be considere a success."
        ),
        type=float,
        default=0.001,
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
