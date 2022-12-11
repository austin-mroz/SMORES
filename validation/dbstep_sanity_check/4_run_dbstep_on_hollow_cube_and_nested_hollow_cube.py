#!python
import argparse
import functools
import pathlib

import dbstep.Dbstep as db
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt

from smores._internal.write_cube import write_cube

matplotlib.use("tkagg")


def main() -> None:
    args = _get_command_line_arguments()
    args.output_directory.mkdir(parents=True, exist_ok=True)

    big_voxel_cube = _get_voxel_cube(5, 16)
    small_voxel_cube = _get_voxel_cube(8, 13)
    joint_voxel_cube = np.logical_or(big_voxel_cube, small_voxel_cube).view(
        np.int8
    )
    if args.plot:
        _plot_voxels(big_voxel_cube)
        _plot_voxels(small_voxel_cube)
        _plot_voxels(joint_voxel_cube)

    write_cube_ = functools.partial(
        write_cube,
        positions=np.array(
            [
                [2.5, 2.5, 2.5],
                [3.2, 3.2, 3.2],
            ],
        ),
        elements=["1", "1"],
        voxel_origin=np.array([-5, -5, -5]),
        voxel_x_vector=np.array([1.0, 0.0, 0.0]),
        voxel_y_vector=np.array([0.0, 1.0, 0.0]),
        voxel_z_vector=np.array([0.0, 0.0, 1.0]),
    )

    big_path = args.output_directory / "big.cube"
    small_path = args.output_directory / "small.cube"
    joint_path = args.output_directory / "joint.cube"

    write_cube_(
        path=big_path,
        voxels=big_voxel_cube,
    )
    write_cube_(
        path=small_path,
        voxels=small_voxel_cube,
    )
    write_cube_(
        path=joint_path,
        voxels=joint_voxel_cube,
    )
    dbstep = functools.partial(
        db.dbstep,
        sterimol=True,
        atom1=0,
        atom2=1,
        surface="density",
        isoval=0.5,
    )
    dbstep(str(big_path))
    dbstep(str(small_path))
    dbstep(str(joint_path))


def _get_voxel_cube(start: int, end: int) -> npt.NDArray:
    voxels = np.zeros((20, 20, 20))
    voxels[start:end, start : end - 1, start] = 1
    voxels[start:end, start : end - 1, end - 1] = 1

    voxels[start:end, start, start : end - 1] = 1
    voxels[start:end, end - 1, start : end - 1] = 1

    voxels[start, start:end, start : end - 1] = 1
    voxels[end - 1, start:end, start : end - 1] = 1
    return voxels


def _plot_voxels(voxels: npt.NDArray) -> None:
    figure = plt.figure()
    axes = figure.add_subplot(projection="3d")
    xs, ys, zs = np.nonzero(voxels)
    axes.scatter(xs, ys, zs)
    plt.show()


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run dbstep on a large voxel cube, a small voxel cube and a "
            "small cube nested inside a large cube."
        ),
    )
    parser.add_argument(
        "-o",
        "--output_directory",
        help="The directory holding the output files.",
        default=pathlib.Path.cwd() / "4_output",
        type=pathlib.Path,
    )
    parser.add_argument(
        "-p",
        "--plot",
        help="If set, a plot of the generated cubes will be displayed.",
        action="store_true",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
