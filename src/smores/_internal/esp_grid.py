from dataclasses import dataclass

import numpy.typing as npt


@dataclass(frozen=True, slots=True)
class ElectrostaticPotentialGrid:
    """
    A 3-D grid of voxels holding electrostatic potentials.

    Attributes:
        grid: The voxel grid.
        voxl_size: The length of a single voxel in each dimension.

    """

    grid: npt.ArrayLike
    voxel_size: tuple[float, float, float]
