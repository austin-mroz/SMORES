import math
from dataclasses import dataclass
from operator import add

import numpy as np
import numpy.typing as npt
import pyvista as pv
import scipy.spatial.distance as distance

from smores._internal.steric_parameters import StericParameters


@dataclass(frozen=True, slots=True)
class BValues:
    Bmin: float
    Bmax: float


def calculate_steric_parameters_from_esp(
    streusel_surface: npt.NDArray[np.float64],
    resolution: np.float32,
    attached_atom_idx: npt.NDArray[np.int_],
    dummy_atom_idx: npt.NDArray[np.int_],
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

    return StericParameters(
        L=l_value,
        B1=b_values.Bmin,
        B5=b_values.Bmax,
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
    return np.average(
        resolution * math.dist(attached_atom_idx, streusel_surface_L_point)
    )


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

    clipped = streusel_surface_point_cloud.clip(
        normal=dummy_atom_idx, origin=attached_atom_idx
    )
    # project clipped surface to plane
    # projected = clipped.project_points_to_plane()

    clipped_and_dummy = pv.PolyData(clipped.points) + pv.PolyData(
        dummy_atom_idx
    )
    merged_projected = clipped_and_dummy.project_points_to_plane(
        origin=attached_atom_idx, normal=dummy_atom_idx
    )

    projected_dummy_atom_idx = merged_projected.points[-1]
    substituent_shadow = merged_projected.points[:-1]

    distances = []
    for point in substituent_shadow:
        distances.append(distance.euclidean(projected_dummy_atom_idx, point))

    b5 = np.max(distances) * resolution**3
    print(f"Bmax: {b5}")  # [0].sum()}")
    index_max = np.argmax(distances)
    max_point = substituent_shadow[index_max]
    circle = pv.Cylinder(
        center=merged_projected.points[-1],
        direction=clipped_and_dummy.points[-1],
        radius=np.max(distances),
        resolution=100,
    )
    circle_arc = pv.CircularArcFromNormal(
        center=merged_projected.points[-1],
        normal=clipped_and_dummy.points[-1],
        polar=substituent_shadow[index_max],
    )

    b1 = _calculate_B1(
        distances, merged_projected.points, dummy_atom_idx, resolution
    )
    return BValues(
        Bmin=b1,
        Bmax=np.average(b5),
    )


def _calculate_B1(
    distances, molecule_shadow, dummy_atom_idx, resolution
) -> np.float64:
    def rotateVector3D(v, theta, axis):
        """Takes a three-dimensional vector v and rotates it by the angle theta around the specified axis."""
        return np.dot(rotationMatrix3D(theta, axis), v)

    def rotationMatrix3D(theta, axis):
        """Return the rotation matrix associated with counterclockwise rotation about
        the given axis by theta radians.
        """
        axis = np.asarray(axis) / np.sqrt(np.dot(axis, axis))
        a = np.cos(theta / 2.0)
        b, c, d = -axis * np.sin(theta / 2.0)
        aa, bb, cc, dd = a**2, b**2, c**2, d**2
        bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
        return np.array(
            [
                [aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc],
            ]
        )

    def drawObject(ax, pts, color="red"):
        """Draws an object on a specified 3D axis with points and lines between consecutive points."""
        map(lambda pt: ax.scatter(*pt, s=10, color=color), pts)
        for k in range(len(pts) - 1):
            x, y, z = zip(*pts[k : k + 2])
            ax.plot(x, y, z, color=color, linewidth=1.0)
        x, y, z = zip(*[pts[-1], pts[0]])
        ax.plot(x, y, z, color=color, linewidth=1.0)

    def normalVector(obj):
        """Takes a set of points, assumed to be flat, and returns a normal vector with unit length."""
        n = np.cross(
            np.array(obj[1]) - np.array(obj[0]),
            np.array(obj[2]) - np.array(obj[0]),
        )
        return n / np.sqrt(np.dot(n, n))

    print("===== CALCULATING B1 =====")
    radius = np.max(distances)
    center_point = molecule_shadow[-1]

    circle_points = get_circle_points(radius, center_point)

    # Set the original object (can be any set of points)
    obj = list(map(tuple, circle_points))
    mObj = tuple(center_point)
    nVecObj = normalVector(obj)

    # Given vector -- molecule_shadow_normal vector
    vec = tuple(center_point - dummy_atom_idx)

    # Find rotation axis and angle.
    rotAxis = normalVector([(0, 0, 0), nVecObj, vec])
    angle = np.arccos(
        np.dot(nVecObj, vec)
        / (np.sqrt(np.dot(vec, vec)) * np.sqrt(np.dot(nVecObj, nVecObj)))
    )

    # Generate the rotated object.
    rotObj = map(lambda pt: rotateVector3D(pt, angle, rotAxis), obj)
    mRotObj = rotateVector3D(mObj, angle, rotAxis)
    rotObj = np.array(list(rotObj))

    nVecRotObj = normalVector(rotObj)

    xadd = center_point[0] - mRotObj[0]
    yadd = center_point[1] - mRotObj[1]
    zadd = center_point[2] - mRotObj[2]

    trotObj = []

    for point in rotObj:
        trotObj.append(list(map(add, point, [xadd, yadd, zadd])))

    # Check if the given vector and the normal of the rotated object are parallel (cross product should be zero).
    b1_vecs = []
    for circle_point in trotObj:
        b1_vecs.append(
            _get_furthest_point_along_vector(
                center_point,
                np.asarray(circle_point),
                molecule_shadow,
                resolution,
            )
        )
    return np.min(b1_vecs)


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
    return resolution**3 * math.dist(starting_point, b1)


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
