#!python

import argparse
import pathlib
from dataclasses import dataclass

import ase.io.cube
import numpy as np
import numpy.typing as npt
import scipy.ndimage

from smores._internal.write_cube import write_cube


def main() -> None:
    args = _get_command_line_arguments()
    args.output_file.parent.mkdir(parents=True, exist_ok=True)

    voxels, atoms = ase.io.cube.read_cube_data(str(args.input_file))
    voxel_grid_params = _get_voxel_grid_params(args.input_file)

    write_cube(
        path=args.output_file,
        voxels=_get_surface(voxels),
        positions=atoms.get_positions() / ase.units.Bohr,
        elements=atoms.get_atomic_numbers(),
        voxel_origin=voxel_grid_params.origin,
        voxel_x_vector=voxel_grid_params.voxel_x_vector,
        voxel_y_vector=voxel_grid_params.voxel_y_vector,
        voxel_z_vector=voxel_grid_params.voxel_z_vector,
    )


@dataclass(frozen=True, slots=True)
class VoxelGridParams:
    origin: npt.NDArray
    voxel_x_vector: npt.NDArray
    voxel_y_vector: npt.NDArray
    voxel_z_vector: npt.NDArray


def _get_voxel_grid_params(cube_file: pathlib.Path) -> VoxelGridParams:
    with open(cube_file) as f:
        next(f)
        next(f)
        origin = np.array(list(map(float, next(f).split()[1:])))
        voxel_x_vector = np.array([float(n) for n in next(f).split()[1:]])
        voxel_y_vector = np.array([float(n) for n in next(f).split()[1:]])
        voxel_z_vector = np.array([float(n) for n in next(f).split()[1:]])

    return VoxelGridParams(
        origin=origin,
        voxel_x_vector=voxel_x_vector,
        voxel_y_vector=voxel_y_vector,
        voxel_z_vector=voxel_z_vector,
    )


def _get_surface(voxels: npt.NDArray[np.float32]) -> npt.NDArray[np.float32]:
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
