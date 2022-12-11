import pathlib

import numpy as np

import smores
from smores._internal.read_cube import read_cube
from smores._internal.write_cube import write_cube


def test_reading_and_writing_is_self_consistent(
    tmp_path: pathlib.Path,
) -> None:
    elements = [1, 2, 3, 4]
    positions = np.array(
        [
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
            [10.0, 11.0, 12.0],
        ]
    )

    rng = np.random.default_rng(4)
    voxel_grid = smores.VoxelGrid(
        voxels=rng.random((15, 15, 15)),
        voxel_origin=rng.random(3),
        voxel_x_vector=rng.random(3),
        voxel_y_vector=rng.random(3),
        voxel_z_vector=rng.random(3),
    )

    cube_path = tmp_path / "molecule.cube"

    write_cube(
        path=cube_path,
        voxels=voxel_grid.voxels,
        positions=positions,
        elements=elements,
        voxel_origin=voxel_grid.voxel_origin,
        voxel_x_vector=voxel_grid.voxel_x_vector,
        voxel_y_vector=voxel_grid.voxel_y_vector,
        voxel_z_vector=voxel_grid.voxel_z_vector,
    )

    cube_data = read_cube(cube_path)

    assert cube_data.atomic_numbers == elements
    assert np.all(np.isclose(cube_data.positions, positions))
    assert np.all(np.isclose(cube_data.voxel_grid.voxels, voxel_grid.voxels))
    assert np.all(
        np.isclose(
            cube_data.voxel_grid.voxel_origin,
            voxel_grid.voxel_origin,
        ),
    )
    assert np.all(
        np.isclose(
            cube_data.voxel_grid.voxel_x_vector,
            voxel_grid.voxel_x_vector,
        ),
    )
    assert np.all(
        np.isclose(
            cube_data.voxel_grid.voxel_y_vector,
            voxel_grid.voxel_y_vector,
        ),
    )
    assert np.all(
        np.isclose(
            cube_data.voxel_grid.voxel_z_vector,
            voxel_grid.voxel_z_vector,
        ),
    )
