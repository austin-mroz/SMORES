import pathlib
from dataclasses import dataclass

import numpy as np
import rdkit.Chem as rdkit

from smores._internal.voxel_grid import VoxelGrid


@dataclass(frozen=True, slots=True)
class CubeFileData:
    atomic_numbers: list[int]
    positions: list[tuple[float, float, float]]
    voxel_grid: VoxelGrid


def read_cube(path: pathlib.Path) -> CubeFileData:
    atomic_numbers: list[int] = []
    positions: list[tuple[float, float, float]] = []
    with open(path) as cube_file:
        _ = next(cube_file)
        _ = next(cube_file)
        num_atoms_, *voxel_origin_ = next(cube_file).split()
        num_atoms = int(num_atoms_)
        voxel_origin = np.array([float(coord) for coord in voxel_origin_])

        num_voxels_x_, *vector_x_ = next(cube_file).split()
        num_voxels_x = int(num_voxels_x_)
        vector_x = np.array([float(coord) for coord in vector_x_])

        num_voxels_y_, *vector_y_ = next(cube_file).split()
        num_voxels_y = int(num_voxels_y_)
        vector_y = np.array([float(coord) for coord in vector_y_])

        num_voxels_z_, *vector_z_ = next(cube_file).split()
        num_voxels_z = int(num_voxels_z_)
        vector_z = np.array([float(coord) for coord in vector_z_])

        for _ in range(num_atoms):
            atom, _, x, y, z = next(cube_file).split()
            atomic_numbers.append(rdkit.Atom(atom).GetAtomicNum())
            positions.append((float(x), float(y), float(z)))

        voxels: list[float] = []
        for line in cube_file.readlines():
            voxels.extend(map(float, line.split()))

    return CubeFileData(
        atomic_numbers=atomic_numbers,
        positions=positions,
        voxel_grid=VoxelGrid(
            voxels=np.array(voxels).reshape(
                num_voxels_x,
                num_voxels_y,
                num_voxels_z,
            ),
            voxel_origin=voxel_origin,
            voxel_x_vector=vector_x,
            voxel_y_vector=vector_y,
            voxel_z_vector=vector_z,
        ),
    )
