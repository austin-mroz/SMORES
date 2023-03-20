import math
from dataclasses import dataclass
from operator import add
import numpy as np
import numpy.typing as npt
import pyvista as pv
import scipy.spatial.distance

from smores._internal.steric_parameters import StericParameters


@dataclass(frozen=True, slots=True)
class BValues:
    Bmin: float
    Bmax: float


def calculate_steric_parameters_from_esp(
        streusel_surface: npt.NDArrayp[np.float64],
        resolution: np.float32,
        attached_atom_idx: npt.NDArray[np.int],
        dummy_atom_idx: npt.NDArray[np.int],
        ) -> StericParameters:
    l_value = _calculate_L(
            attached_atom_idx,
            dummy_atom_idx,
            streusel_surface,
            resolution,
            )
    b_values = _calculate_B(attached_atom_idx, dummy_atom_idx, streusel_surface, resolution,)

    return StericParameters(
            L=l_value,
            B1=b_values.Bmin,
            B5=b_values.Bmax
            )


def _calculate_L(
    attached_atom_idx: npt.NDArray,  # substituent
    dummy_atom_idx: npt.NDArray,  # core
    streusel_surface: npt.NDArray[np.float64],
    resolution: np.float32,
) -> np.float32:
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

    return resolution * math.dist(attached_atom_idx, streusel_surface_L_point)

def _calculate_B(
    attached_atom_idx: npt.NDArray,
    dummy_atom_idx: npt.NDArray,
    streusel_surface: npt.NDArray[np.float64],
    resolution: float,
) -> BValues:
    streusel_surface_idx = np.argwhere(streusel_surface)

    streusel_surface_point_cloud = pv.PolyData(streusel_surface_idx)

    # define plane with center at core and normal along substituent
    b_vector_plane = pv.Plane(attached_atom_idx, dummy_atom_idx)

    clipped = streusel_surface_point_cloud.clip(normal=dummy_atom_idx, origin=attached_atom_idx)

    # project clipped surface to plane
    projected = clipped.project_points_to_plane()
    projected.plot()

    clipped_and_dummy = pv.PolyData(clipped.points) + pv.PolyData(dummy_atom_idx)
    merged_projected = clipped_and_dummy.project_points_to_plane(origin=attached_atom_idx, normal=dummy_atom_idx)
    p = pv.Plotter()
    p.add_mesh(clipped, color="w")
    p.add_mesh(pv.PolyData(dummy_atom_idx), color="b", point_size=20)
    p.add_mesh(clipped_and_dummy, color="g")
    p.add_mesh(pv.PolyData(clipped_and_dummy.points[-1]), point_size=30, color="#69FAAB")
    p.add_mesh(merged_projected)
    p.add_mesh(pv.PolyData(merged_projected.points[-1]), point_size=20, color='#69FAAB')
    p.show()

    projected_dummy_atom_idx = merged_projected.points[-1]
    substituent_shadow = merged_projected.points[:-1]

    distances = []
    for point in substituent_shadow:
        distances.append(distance.euclidean(projected_dummy_atom_idx, point))

    b5 = np.max(distances)*resolution**3
    print(f"Bmax: {b5}") # [0].sum()}")
    index_max = np.argmax(distances)
    max_point = substituent_shadow[index_max]
    circle = pv.Cylinder(center=merged_projected.points[-1], direction=clipped_and_dummy.points[-1], radius=np.max(distances), resolution=100)
    circle_arc = pv.CircularArcFromNormal(center=merged_projected.points[-1],normal=clipped_and_dummy.points[-1], polar=substituent_shadow[index_max])

    b1 = calculate_B1(distances, merged_projected.points, dummy_atom_idx, resolution)
    print(f"Bmin: {b1}")

    return BValues(Bmin=b1, Bmax=b5)



def _get_furthest_point_along_vector(
        starting_point: npt.NDArray,
        ending_point: npt.NDArray,
        molecule_shadow: npt.NDArray,
        ) -> np.float64:
    product = np.cross(molecule_shadow - ending_point,
                       starting_point - ending_point,
                       )
    if product.ndim == 2:
        distances = np.linalg.norm(product, axis=1)
    else:
        distances = np.abs(product)
    b1 = molecule_shadow[distances.argmax()]
    return 0.2 * math.dist(starting_point, b1)


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

