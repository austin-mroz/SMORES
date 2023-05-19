#!python
import argparse
import pathlib

import atomlite
import flour
import morfeus
import numpy as np
import numpy.typing as npt
import pyvista as pv
import scipy

import smores


def main() -> None:
    args = _get_command_line_arguments()

    radii_types = ("alvarez", "bondi", "crc", "rahm", "pyykko", "truhlar")

    database = atomlite.Database(args.database)

    for system in database.get_entries():
        print(system.key)

        smores_molecule = smores.Molecule.from_xyz_file(
            path=pathlib.Path(system.properties["xyz_file"]),
            dummy_index=system.properties["dummy_index"],
            attached_index=system.properties["attached_index"],
        )
        smores_params = smores_molecule.get_steric_parameters()
        new_entry = atomlite.PropertyEntry(
            key=system.key,
            properties={
                "streusel_L": smores_params.L,
                "streusel_B1": smores_params.B1,
                "streusel_B5": smores_params.B5,
            },
        )
        database.update_properties(new_entry)

        for radii_type in radii_types:
            sterimol = _get_sterimol(system, radii_type)
            new_entry = atomlite.PropertyEntry(
                key=system.key,
                properties={
                    f"{radii_type}_L": sterimol.L_value,
                    f"{radii_type}_B1": sterimol.B_1_value,
                    f"{radii_type}_B5": sterimol.B_5_value,
                },
            )
            database.update_properties(new_entry)

        smores_esp_molecule = smores.EspMolecule.from_cube_file(
            path=pathlib.Path(system.properties["esp_file"]),
            dummy_index=int(system.properties["dummy_index"]),
            attached_index=int(system.properties["attached_index"]),
        )

        esp_smores_params = smores_esp_molecule.get_steric_parameters()

        new_entry = atomlite.PropertyEntry(
            key=system.key,
            properties={
                "streusel_cube_L": esp_smores_params.L * 0.529,
                "streusel_cube_B1": esp_smores_params.B1 * 0.529,
                "streusel_cube_B5": esp_smores_params.B5 * 0.529,
            },
        )
        database.update_properties(new_entry)
    database.connection.commit()


def _get_sterimol(
    molecule: atomlite.Entry, radii_type: str
) -> morfeus.Sterimol:
    elements, coordinates = morfeus.read_xyz(molecule.properties["xyz_file"])
    return morfeus.Sterimol(
        elements=elements,
        coordinates=coordinates,
        dummy_index=molecule.properties["dummy_index"] + 1,
        attached_index=molecule.properties["attached_index"] + 1,
        radii_type=radii_type,
    )


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Calculate and SMORES and sterimol parameters.",
    )

    parser.add_argument(
        "-d",
        "--database",
        help=(
            'An atomlite database file with properties: "core", "substituent", '
            '"xyz_file", "esp_file", dummy_index" and "attached_index".'
        ),
        type=pathlib.Path,
        default=pathlib.Path.cwd() / "catalyst_systems.db",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
