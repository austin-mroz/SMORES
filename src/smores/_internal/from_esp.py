import math
import pathlib
from dataclasses import dataclass
from operator import add

import numpy as np
import numpy.typing as npt
import pyvista as pv
import scipy.spatial.distance as distance

from smores._internal.steric_parameters import StericParameters


@dataclass(frozen=True, slots=True)
class Lvalue:
    L: float
    L_idx: npt.NDArray[np.float64]


@dataclass(frozen=True, slots=True)
class BValues:
    Bmin: float
    Bmax: float
    Bmin_idx: npt.NDArray[np.float64]
    Bmax_idx: npt.NDArray[np.float64]


def calculate_steric_parameters_from_esp(
    streusel_surface: npt.NDArray[np.float64],
    resolution: np.float32,
    attached_atom_idx: npt.NDArray[np.int_],
    dummy_atom_idx: npt.NDArray[np.int_],
    plot: bool = False,
    output_path: None | pathlib.Path = None,
) -> StericParameters:
    l_value = _calculate_L(
        attached_atom_idx,
        dummy_atom_idx,
        streusel_surface,
        resolution.sum(axis=1),
    )

    b_values = _calculate_B(
        attached_atom_idx,
        dummy_atom_idx,
        streusel_surface,
        resolution.sum(axis=1),
    )

    if plot:
        plot_steric_parameters(
            attached_atom_idx,
            dummy_atom_idx,
            streusel_surface,
            l_value,
            b_values,
            output_path,
        )

    return StericParameters(
        L=l_value.L,
        B1=b_values.Bmin,
        B5=b_values.Bmax,
    )


def plot_steric_parameters(
    attached_atom_idx: npt.NDArray,
    dummy_atom_idx: npt.NDArray,
    streusel_surface: npt.NDArray,
    l_value: Lvalue,
    b_values: BValues,
    output_path: pathlib.Path,
) -> None:
    pink = "#ff79c6"
    green = "#50fa7b"
    cyan = "#8be9fd"
    orange = "#ffb86c"
    purple = "#bd93f9"
    dustyblue = "#6277a4"

    streusel_surface_idx = np.argwhere(streusel_surface)

    clipped = pv.PolyData(streusel_surface_idx).clip(
        normal=attached_atom_idx - dummy_atom_idx,
        origin=attached_atom_idx,
    )

    shadow = clipped.project_points_to_plane(
        normal=attached_atom_idx - dummy_atom_idx,
        origin=attached_atom_idx,
    )

    p = pv.Plotter(off_screen=True, shape=(1, 2))

    p.set_background("white")

    p.add_mesh(
        pv.PolyData(attached_atom_idx),
        point_size=20,
        color=pink,
    )
    p.add_mesh(
        pv.PolyData(dummy_atom_idx),
        point_size=20,
        color=cyan,
    )
    p.add_mesh(
        clipped,
        opacity=0.3,
    )
    p.add_mesh(
        shadow,
        color=purple,
        opacity=0.3,
    )

    p.add_mesh(
        pv.PolyData(b_values.Bmin_idx),
        point_size=20,
        color=green,
    )
    p.add_mesh(
        pv.PolyData(b_values.Bmax_idx),
        point_size=20,
        color=green,
    )
    p.add_mesh(
        pv.PolyData(l_value.L_idx),
        point_size=20,
        color=green,
    )

    p.view_vector(
        vector=attached_atom_idx - dummy_atom_idx,
        viewup=b_values.Bmax_idx - b_values.Bmin_idx,
    )

    p.subplot(0, 1)
    p.set_background("white")
    p.add_mesh(
        pv.PolyData(attached_atom_idx),
        point_size=20,
        color=pink,
    )
    p.add_mesh(
        pv.PolyData(dummy_atom_idx),
        point_size=20,
        color=cyan,
    )
    p.add_mesh(
        clipped,
        opacity=0.3,
        color=dustyblue,
    )
    p.add_mesh(
        shadow,
        color=purple,
        opacity=0.3,
    )

    p.add_mesh(
        pv.PolyData(b_values.Bmin_idx),
        point_size=20,
        color=green,
    )
    p.add_mesh(
        pv.PolyData(b_values.Bmax_idx),
        point_size=20,
        color=green,
    )
    p.add_mesh(
        pv.PolyData(l_value.L_idx),
        point_size=20,
        color=green,
    )

    p.view_vector(
        vector=b_values.Bmax_idx - b_values.Bmin_idx,
        viewup=attached_atom_idx - dummy_atom_idx,
    )

    p.screenshot(output_path)


def _calculate_L(
    attached_atom_idx: npt.NDArray,  # substituent
    dummy_atom_idx: npt.NDArray,  # core
    streusel_surface: npt.NDArray[np.float64],
    resolution: np.float32,
) -> Lvalue:
    streusel_surface_idx = np.argwhere(streusel_surface)

    clipped = pv.PolyData(streusel_surface_idx).clip(
        normal=attached_atom_idx - dummy_atom_idx,
        origin=attached_atom_idx,
    )

    product = np.cross(
        np.array(clipped.points) - dummy_atom_idx,
        attached_atom_idx - dummy_atom_idx,
    )
    if product.ndim == 2:
        distances = np.linalg.norm(product, axis=1)
    else:
        distances = np.abs(product)

    streusel_surface_L_point = clipped.points[distances.argmin()]
    return Lvalue(
        L=np.average(
            resolution * math.dist(attached_atom_idx, streusel_surface_L_point)
        ),
        L_idx=np.array(streusel_surface_L_point),
    )


def _calculate_B(
    attached_atom_idx: npt.NDArray,
    dummy_atom_idx: npt.NDArray,
    streusel_surface: npt.NDArray[np.float64],
    resolution: float,
) -> BValues:
    streusel_surface_idx = np.argwhere(streusel_surface)

    clipped = pv.PolyData(streusel_surface_idx).clip(
        normal=attached_atom_idx - dummy_atom_idx,
        origin=attached_atom_idx,
    )

    shadow = clipped.project_points_to_plane(
        normal=attached_atom_idx - dummy_atom_idx,
        origin=attached_atom_idx,
    )

    shadow_surface = shadow.delaunay_2d(alpha=1.0)
    shadow_edges = shadow_surface.extract_feature_edges(
        feature_angle=30,
        boundary_edges=True,
        non_manifold_edges=False,
        feature_edges=False,
        manifold_edges=False,
        clear_data=True,
    )

    b_value_distances = []
    for point in shadow_edges.points:
        b_value_distances.append(
            distance.euclidean(
                point,
                attached_atom_idx,
            )
        )
    b1_idx = shadow_edges.points[np.array(b_value_distances).argmin()]
    b5_idx = shadow_edges.points[np.array(b_value_distances).argmax()]

    b1 = min(b_value_distances)
    b5 = max(b_value_distances)

    return BValues(
        Bmin=b1 * np.average(resolution),
        Bmax=b5 * np.average(resolution),
        Bmin_idx=np.array(b1_idx),
        Bmax_idx=np.array(b5_idx),
    )


def _get_furthest_point_along_vector(
    starting_point: npt.NDArray,
    ending_point: npt.NDArray,
    molecule_shadow: npt.NDArray,
    resolution: np.float32,
) -> np.float64:
    product = np.cross(
        molecule_shadow - ending_point,
        starting_point - ending_point,
    )
    if product.ndim == 2:
        distances = np.linalg.norm(product, axis=1)
    else:
        distances = np.abs(product)
    b1 = molecule_shadow[distances.argmax()]
    return resolution * math.dist(starting_point, b1)


def get_circle_points(radius, center_point):
    return np.asarray(
        [
            [
                center_point[0] + radius * math.cos(2 * math.pi / 1000 * x),
                center_point[1] + radius * math.sin(2 * math.pi / 1000 * x),
                center_point[2],
            ]
            for x in range(0, 1000)
        ]
    )


def calculate_angle_between(circle, molecule) -> float:
    circle_unit_vector = circle / np.linalg.norm(circle)
    molecule_unit_vector = molecule / np.linalg.norm(molecule)

    return np.arccos(
        np.clip(np.dot(circle_unit_vector, molecule_unit_vector), -1.0, 1.0)
    )


def calc_distance(projected_dummy_atom_idx, point):
    return distance.euclidean(projected_dummy_atom_idx, point)


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
