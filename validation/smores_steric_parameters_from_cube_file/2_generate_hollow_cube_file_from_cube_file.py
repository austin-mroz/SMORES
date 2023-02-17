#!python

import argparse
import pathlib

import flour
import numpy as np
import numpy.typing as npt
import scipy.ndimage


def main() -> None:
    args = _get_command_line_arguments()
    args.output_file.parent.mkdir(parents=True, exist_ok=True)

    cube_data = flour.read_cube(args.input_file)

    flour.write_cube(
        path=args.output_file,
        title1="t1",
        title2="t2",
        charges=cube_data.charges,
        voxels=_get_surface(cube_data.grid.voxels).astype(float),
        positions=cube_data.positions,
        atoms=cube_data.atoms,
        voxel_origin=cube_data.grid.origin,
        voxel_size=cube_data.grid.voxel_size,
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
