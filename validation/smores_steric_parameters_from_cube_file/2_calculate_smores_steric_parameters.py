#!python

import argparse
import pathlib

import flour
import numpy as np
import numpy.typing as npt
import pyvista as pv
import scipy.ndimage


def main() -> None:
    args = _get_command_line_arguments()
    args.output_file.parent.mkdir(parents=True, exist_ok=True)

    cube_data = flour.read_cube(args.input_file)

    print(type(cube_data))
    cube_data_positions_idx = _convert_euclidean_positions_to_indices(
        cube_data
    )

    streusel_surface = _get_surface(cube_data.grid.voxels)

    print(streusel_surface.shape)
    _calculate_L(
        cube_data_positions_idx[0],
        cube_data_positions_idx[1],
        streusel_surface,
    )


def _calculate_L(
    attached_atom_idx: npt.NDArray,  # substituent
    dummy_atom_idx: npt.NDArray,  # core
    streusel_surface: npt.NDArray[np.float64],
    # cube_data_positions_idx: npt.NDArray[np.float64],
) -> float:
    # vector along which L terminal point should lie on the surface
    vector_L = attached_atom_idx - dummy_atom_idx

    # extract surface points
    streusel_surface_idx = np.argwhere(streusel_surface)

    streusel_surface_point_cloud = pv.PolyData(streusel_surface_idx)
    print(streusel_surface_point_cloud)

    # L vector we care about
    line_L = pv.Line(dummy_atom_idx, attached_atom_idx)

    plane = pv.Plane()

    plane = plane.interpolate(streusel_surface_point_cloud)
    sample = plane.sample_over_line(dummy_atom_idx, attached_atom_idx)
    print(sample)


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
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input_file",
        type=pathlib.Path,
        default=pathlib.Path.cwd() / "1_output" / "ESP.cube",
    )
    parser.add_argument(
        "-o",
        "--output_file",
        type=pathlib.Path,
        default=pathlib.Path.cwd() / "2_output" / "hollow_ESP.cube",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
