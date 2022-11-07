from dataclasses import dataclass

import numpy.typing as npt


@dataclass(frozen=True, slots=True)
class VoxelGrid:
    """
    A 3-D grid of voxels holding electrostatic potentials.

    """

    #: The voxels of the grid, represented as a 3-D array.
    voxels: npt.NDArray
    #: The length of a single voxel along the x, y and z dimensions.
    voxel_size: npt.NDArray
    #: The origin of voxels.
    voxel_origin: npt.NDArray
