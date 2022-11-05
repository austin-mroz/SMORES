from dataclasses import dataclass

import numpy as np
import numpy.typing as npt


@dataclass(frozen=True, slots=True)
class ElectrostaticPotentialGrid:
    """
    A 3-D grid of voxels holding electrostatic potentials.

    """

    #: The voxel grid, represented as a 3-D array.
    grid: npt.NDArray[np.float32]
    #: The length of a single voxel along the x, y and z dimensions.
    voxel_size: tuple[float, float, float]
