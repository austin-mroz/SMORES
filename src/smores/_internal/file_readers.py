import pathlib
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt


@dataclass(frozen=True, slots=True)
class XyzData:
    atoms: list[str]
    positions: npt.NDArray[np.float32]


def read_xyz(
    path: pathlib.Path,
) -> XyzData:

    with open(path, "r") as f:
        lines = f.readlines()

    elements = []
    coordinates = []
    for line in lines[2:]:
        element, x, y, z = line.strip().split()
        elements.append(element)
        coordinates.append((float(x), float(y), float(z)))

    return XyzData(
        elements=elements,
        coordinates=np.array(coordinates),
    )
