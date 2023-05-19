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

        print(system.key)
        smores_esp_molecule = smores.EspMolecule.from_cube_file(
            path=pathlib.Path(system.properties["esp_file"]),
            dummy_index=int(system.properties["dummy_index"]),
            attached_index=int(system.properties["attached_index"]),
        )

        esp_smores_params = smores_esp_molecule.get_steric_parameters()

        new_entry = atomlite.PropertyEntry(
            key=system.key,
            properties={
                "streusel_cube_L": esp_smores_params.L,
                "streusel_cube_B1": esp_smores_params.B1,
                "streusel_cube_B5": esp_smores_params.B5,
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


def _plot_surface(
    cube_file: pathlib.Path,
    attached_atom_idx: int,
    dummy_atom_idx: int,
) -> None:
    cube_data = flour.read_cube(cube_file)

    cube_data_positions_idx = _convert_euclidean_positions_to_indices(
        cube_data
    )

    streusel_surface = _get_surface(cube_data.grid.voxels)

    streusel_surface_idx = np.argwhere(streusel_surface)
    streusel_surface_point_cloud = pv.PolyData(streusel_surface_idx)

    p = pv.Plotter()
    p.add_mesh(streusel_surface_point_cloud)
    p.add_mesh(
        pv.PolyData(cube_data_positions_idx[dummy_atom_idx]),
        point_size=20,
        color="#69FAAB",
    )
    p.add_mesh(
        pv.PolyData(cube_data_positions_idx[attached_atom_idx]),
        point_size=20,
        color="#69FAAB",
    )
    p.show(screenshot="c_h.png")


def _convert_euclidean_positions_to_indices(
    cube_data,  # : flour.CubeData,
) -> npt.NDArray:
    # translate origin and positions to 0,0,0
    translated_position_matrix = cube_data.positions - cube_data.grid.origin

    lattice_lengths = (
        cube_data.grid.voxel_size.sum(axis=0) * cube_data.grid.voxels.shape
    )

    percent_along_lattice_vector = np.divide(
        translated_position_matrix, lattice_lengths
    )

    return (percent_along_lattice_vector * cube_data.grid.voxels.shape).astype(
        int
    )


def _get_surface(voxels: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    cv = 0.0004
    is_vacuum = voxels < cv
    is_non_vacuum = np.logical_not(is_vacuum)

    weights = np.zeros((3, 3, 3))
    weights[1, 1, 2] = 1
    weights[1, 1, 0] = 1
    weights[:, :, 1] = [
        [0, 1, 0],
        [1, 0, 1],
        [0, 1, 0],
    ]

    convolution_result = scipy.ndimage.convolve(
        input=is_vacuum.view(np.int8),
        weights=weights,
        mode="wrap",
    )
    np.multiply(
        convolution_result,
        is_non_vacuum.view(np.int8),
        out=convolution_result,
    )
    return (convolution_result > 0).view(np.int8)


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
        default=pathlib.Path.cwd() / "common_carbon_substituents.db",
    )

    return parser.parse_args()


if __name__ == "__main__":
    main()
