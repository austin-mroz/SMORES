import pathlib
import typing

import numpy.typing as npt


def write_cube(
    path: pathlib.Path,
    voxels: npt.NDArray,
    positions: npt.NDArray,
    elements: typing.Sequence[str],
    voxel_origin: npt.NDArray,
    voxel_dimensions: npt.NDArray,
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
            elemental symbol.

        voxel_origin:
            Origin of the voxels.

        voxel_dimensions:
            The length of a single voxel in the x, y and
            z dimensions.

    """

    origin_x, origin_y, origin_z = voxel_origin
    voxel_x_length, voxel_y_length, voxel_z_length = voxel_dimensions
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
                f"{num_voxels_x: >5} {voxel_x_length: >11.6f} {0.: >11.6f} "
                f"{0.: >11.6f}"
            ),
            (
                f"{num_voxels_y: >5} {0.: >11.6f} {voxel_y_length: >11.6f} "
                f"{0.: >11.6f}"
            ),
            (
                f"{num_voxels_z: >5} {0.: >11.6f} {0.: >11.6f} "
                f"{voxel_z_length: >11.6f}"
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
