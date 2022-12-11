import pathlib
from collections import abc

import numpy.typing as npt


def write_cube(
    path: pathlib.Path,
    voxels: npt.NDArray,
    positions: npt.NDArray,
    elements: abc.Collection[str | int],
    voxel_origin: npt.NDArray,
    voxel_x_vector: npt.NDArray,
    voxel_y_vector: npt.NDArray,
    voxel_z_vector: npt.NDArray,
) -> None:
    """
    Write a ``.cube`` file.

    Paramters:

        path:
            The path to the file being written.

        voxels:
            A 3-D grid contaning voxel values.

        positions:
            An N x 3 matrix of atomic coordinates.

        elements:
            For each element of the molecule, its
            elemental symbol or atomic number.

        voxel_origin:
            Origin of the voxels.

        voxel_x_vector:
            The x vector of a single voxel.

        voxel_y_vector:
            The y vector of a single voxel.

        voxel_z_vector:
            The z vector of a single voxel.

    """

    origin_x, origin_y, origin_z = voxel_origin
    num_voxels_x, num_voxels_y, num_voxels_z = voxels.shape
    with open(path, "w") as cube:
        lines = [
            " title",
            " title2",
            (
                f"{len(elements): >5} {origin_x: >11.6f} {origin_y: >11.6f} "
                f"{origin_z: >11.6f}"
            ),
            (
                f"{num_voxels_x: >5} {voxel_x_vector[0]: >11.6f} "
                f"{voxel_x_vector[1]: >11.6f} "
                f"{voxel_x_vector[2]: >11.6f}"
            ),
            (
                f"{num_voxels_y: >5} {voxel_y_vector[0]: >11.6f}"
                f"{voxel_y_vector[1]: >11.6f} "
                f"{voxel_y_vector[2]: >11.6f}"
            ),
            (
                f"{num_voxels_z: >5} {voxel_z_vector[0]: >11.6f} "
                f"{voxel_z_vector[1]: >11.6f} "
                f"{voxel_z_vector[2]: >11.6f}"
            ),
        ]
        for element, [x, y, z] in zip(elements, positions):
            lines.append(
                f"{element: >5} {0.: >11.6f} {x: >11.6f} "
                f"{y: >11.6f} {z: >11.6f}"
            )

        num_columns = 6
        for i in range(num_voxels_x):
            for j in range(num_voxels_y):
                line = []
                for k in range(num_voxels_z):
                    line.append(f"{voxels[i, j, k]: >12.5E}")
                    if len(line) == num_columns:
                        lines.append(_pad_left(" ".join(line)))
                        line = []
                lines.append(_pad_left(" ".join(line)))
                line = []
        content = "\n".join(lines)
        cube.write(f"{content}\n")


def _pad_left(s: str) -> str:
    return f" {s}"
