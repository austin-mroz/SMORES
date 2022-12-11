from dataclasses import dataclass

import numpy.typing as npt


@dataclass(frozen=True, slots=True)
class VoxelGrid:
    """
    A 3-D grid of voxels holding electrostatic potentials.

    """

    #: The voxels of the grid, represented as a 3-D array.
    voxels: npt.NDArray
    #: The origin of the voxels.
    voxel_origin: npt.NDArray
    #: The x vector of a single voxel.
    voxel_x_vector: npt.NDArray
    #: The y vector of a single voxel.
    voxel_y_vector: npt.NDArray
    #: The z vector of a single voxel.
    voxel_z_vector: npt.NDArray
