#!python

import argparse
import math
import pathlib

import flour
import numpy as np
import numpy.typing as npt
import pyvista as pv
import scipy.ndimage
from scipy.spatial import distance
from sklearn.preprocessing import minmax_scale


def main() -> None:
    args = _get_command_line_arguments()
    args.output_file.parent.mkdir(parents=True, exist_ok=True)

    cube_data = flour.read_cube(args.input_file)

    print(type(cube_data))
    cube_data_positions_idx = _convert_euclidean_positions_to_indices(
        cube_data
    )

    streusel_surface = _get_surface(cube_data.grid.voxels)
    """
    print(streusel_surface.shape)
    _calculate_L(
        cube_data_positions_idx[0],
        cube_data_positions_idx[1],
        streusel_surface,
    )
    _plot_L(cube_data_positions_idx[0],
            cube_data_positions_idx[1],
            streusel_surface,)
    resolution = np.sum(cube_data.grid.voxel_size, axis=1)[0]
    """
    _calculate_B(
            cube_data_positions_idx[0],
            cube_data_positions_idx[1],
            streusel_surface,
            cube_data.grid.voxel_size,
            )


def _calculate_B(
    attached_atom_idx: npt.NDArray,
    dummy_atom_idx: npt.NDArray,
    streusel_surface: npt.NDArray[np.float64],
    resolution: float,
) -> None:
    streusel_surface_idx = np.argwhere(streusel_surface)

    streusel_surface_point_cloud = pv.PolyData(streusel_surface_idx)

    # define plane with center at core and normal along substituent
    b_vector_plane = pv.Plane(attached_atom_idx, dummy_atom_idx)

    clipped = streusel_surface_point_cloud.clip(normal=dummy_atom_idx, origin=attached_atom_idx)

    p = pv.Plotter()
    # p.add_mesh(streusel_surface_point_cloud, color="w")
    p.add_mesh(clipped, color="b")
    #p.show()

    # project clipped surface to plane
    projected = clipped.project_points_to_plane()
    projected.plot()

    print(type(projected))

    print(clipped.points)
    print(projected.points)

    print(dummy_atom_idx)
    print(attached_atom_idx)

    p = pv.Plotter()
    p.add_mesh(clipped)
    p.add_mesh(pv.PolyData(dummy_atom_idx), point_size=20, color='#69FAAB')
    #p.show()

    clipped_and_dummy = pv.PolyData(clipped.points) + pv.PolyData(dummy_atom_idx)
    merged_projected = clipped_and_dummy.project_points_to_plane(origin=attached_atom_idx, normal=dummy_atom_idx)
    p = pv.Plotter()
    p.add_mesh(clipped, color="w")
    p.add_mesh(pv.PolyData(dummy_atom_idx), color="b", point_size=20)
    p.add_mesh(clipped_and_dummy, color="g")
    p.add_mesh(pv.PolyData(clipped_and_dummy.points[-1]), point_size=30, color="#69FAAB")
    p.add_mesh(merged_projected)
    p.add_mesh(pv.PolyData(merged_projected.points[-1]), point_size=20, color='#69FAAB')
    #p.show()

    projected_dummy_atom_idx = merged_projected.points[-1]
    substituent_shadow = merged_projected.points[:-1]

    distances = []
    for point in substituent_shadow:
        distances.append(distance.euclidean(projected_dummy_atom_idx, point))

    b5 = np.max(distances)*resolution
    print(b5[0].sum())
    print(np.max(distances))
    index_max = np.argmax(distances)
    max_point = substituent_shadow[index_max]
    print(index_max)
    print(max_point)
    circle = pv.Cylinder(center=merged_projected.points[-1], direction=clipped_and_dummy.points[-1], radius=np.max(distances), resolution=100)
    circle_arc = pv.CircularArcFromNormal(center=merged_projected.points[-1],normal=clipped_and_dummy.points[-1], polar=substituent_shadow[index_max])

    p = pv.Plotter()
    p.add_mesh(merged_projected)
    p.add_mesh(pv.PolyData(merged_projected.points[-1]), point_size=20, color='#69FAAB')
    # p.add_mesh(circle)
    p.add_mesh(circle_arc, color='#FFC0CB', line_width=10)
    p.show()

    r = np.max(distances)
    circle_circumference_points = [(math.cos(2*pi/1000*x)*r,math.sin(2*pi/1000*x)*r) for x in range(0,1000)]


def calc_distance(projected_dummy_atom_idx, point):
    return distance.euclidean(projected_dummy_atom_idx, point)

def _calculate_L(
    attached_atom_idx: npt.NDArray,  # substituent
    dummy_atom_idx: npt.NDArray,  # core
    streusel_surface: npt.NDArray[np.float64],
) -> np.float64:
    streusel_surface_idx = np.argwhere(streusel_surface)
    product = np.cross(
        streusel_surface_idx - dummy_atom_idx,
        attached_atom_idx - dummy_atom_idx,
    )
    if product.ndim == 2:
        distances = np.linalg.norm(product, axis=1)
    else:
        distances = np.abs(product)

    streusel_surface_L_point = streusel_surface_idx[distances.argmin()]

    return 0.2 * math.dist(attached_atom_idx, streusel_surface_L_point)


def _plot_L(
    attached_atom_idx: npt.NDArray,  # substituent
    dummy_atom_idx: npt.NDArray,  # core
    streusel_surface: npt.NDArray[np.float64],
) -> None:
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

    # normalize vector_L
    u = np.divide(vector_L, np.linalg.norm(vector_L))

    terminal_line_L_point = attached_atom_idx + 100 * u

    line_L = pv.Line(dummy_atom_idx, terminal_line_L_point)
    sample = plane.sample_over_line(dummy_atom_idx, terminal_line_L_point)
    p = pv.Plotter()
    p.add_mesh(streusel_surface_point_cloud, style="wireframe", color="w")
    p.add_mesh(line_L)
    p.add_mesh(plane)
    p.show()

    streusel_surface_point_cloud.compute_implicit_distance(plane, inplace=True)
    dist = streusel_surface_point_cloud["implicit_distance"]
    p = pv.Plotter()
    p.add_mesh(streusel_surface_point_cloud, scalars="implicit_distance")
    p.add_mesh(line_L)
    p.show()

    streusel_surface_point_cloud["values"] = np.full(
        streusel_surface_idx.shape[0], 1
    )
    streusel_surface_point_cloud.plot_over_line(
        dummy_atom_idx, terminal_line_L_point, scalars="values"
    )


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
