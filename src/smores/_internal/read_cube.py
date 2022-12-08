import pathlib
from dataclasses import dataclass

import numpy as np

from smores._internal.voxel_grid import VoxelGrid


@dataclass(frozen=True, slots=True)
class Atom:
    element: int
    position: tuple[float, float, float]


@dataclass(frozen=True, slots=True)
class CubeFileData:
    atoms: list[Atom]
    voxel_grid: VoxelGrid


def read_cube_file(path: pathlib.Path) -> CubeFileData:
    atoms: list[Atom] = []
    with open(path) as cube_file:
        _ = next(cube_file)
        _ = next(cube_file)
        num_atoms, *voxel_origin = next(cube_file).split()
        num_atoms = int(num_atoms)
        voxel_origin = np.array([float(coord) for coord in voxel_origin])

        num_voxels_x, *vector_x = next(cube_file).split()
        num_voxels_x = int(num_voxels_x)
        vector_x = np.array([float(coord) for coord in vector_x])

        num_voxels_y, *vector_y = next(cube_file).split()
        num_voxels_y = int(num_voxels_y)
        vector_y = np.array([float(coord) for coord in vector_y])

        num_voxels_z, *vector_z = next(cube_file).split()
        num_voxels_z = int(num_voxels_z)
        vector_z = np.array([float(coord) for coord in vector_z])

        for i in range(num_atoms):
            atomic_number, _, *position = next(cube_file).split()
            atoms.append(
                Atom(
                    element=int(atomic_number),
                    position=tuple(map(float, position)),
                )
            )
        for i in range(num_voxels_x):
            for j in range(num_voxels_y):
                for k in range(num_voxels_z):
                    pass
