import argparse

import dbstep.Dbstep as db


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("log_file")
    parser.add_argument("cube_file")
    return parser.parse_args()


def main() -> None:

    cli_args = _get_command_line_arguments()
    classic_mol = db.dbstep(
        cli_args.log_file,
        atom1=2,
        atom2=5,
        commandline=True,
        verbose=True,
        sterimol=True,
        measure="classic",
    )

    print(
        f"CLASSIC "
        f"-- L: {classic_mol.L} "
        f"-- BMIN: {classic_mol.Bmin} "
        f"-- BMAX: {classic_mol.Bmax}"
    )

    grid_mol = db.dbstep(
        cli_args.cube_file,
        atom1=2,
        atom2=5,
        commandline=True,
        verbose=True,
        sterimol=True,
        measure="grid",
        surface="vdw",
    )

    print(
        f"GRID "
        f"-- L: {grid_mol.L} "
        f"-- BMIN: {grid_mol.Bmin} "
        f"-- BMAX: {grid_mol.Bmax}"
    )


if __name__ == "__main__":
    main()
