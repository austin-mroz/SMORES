#!python

import argparse
import pathlib

import dbstep.Dbstep as db


def main() -> None:
    args = _get_command_line_arguments()

    db.dbstep(
        str(args.cube_file),
        sterimol=True,
        atom1=2,
        atom2=5,
        surface="density",
        isoval=0.0004,
    )
    db.dbstep(
        str(args.hollow_cube_file),
        sterimol=True,
        atom1=2,
        atom2=5,
        surface="density",
    )


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run dbstep on a regular cube file as well as one where only "
            "the surface voxels are non-zero."
        ),
    )
    parser.add_argument(
        "--cube_file",
        type=pathlib.Path,
        default=pathlib.Path.cwd() / "1_output" / "ESP.cube",
    )
    parser.add_argument(
        "--hollow_cube_file",
        type=pathlib.Path,
        default=pathlib.Path.cwd() / "2_output" / "hollow_ESP.cube",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
