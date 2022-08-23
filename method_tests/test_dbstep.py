import argparse

import dbstep.Dbstep as db


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("log_file")
    parser.add_argument("cube_file")
    return parser.parse_args()


def main() -> None:

    cli_args = _get_command_line_arguments()
    db.dbstep(
        cli_args.log_file,
        atom1=2,
        atom2=5,
        commandline=True,
        verbose=True,
        sterimol=True,
        measure="classic",
    )

    db.dbstep(
        cli_args.cube_file,
        atom1=2,
        atom2=5,
        commandline=True,
        verbose=True,
        sterimol=True,
        measure="grid",
        surface="vdw",
    )


if __name__ == "__main__":
    main()
