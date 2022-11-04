import pathlib
import typing
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

from . import molecule
from . import constants


@dataclass(frozen=True)
class XyzData:
    elements: list[str]
    coordinates: npt.NDArray[np.float32]


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
        coordinates.append((x, y, z))
    return XyzData(
        elements=elements,
        coordinates=np.array(coordinates).astype(np.float32),
    )


def write_xyz(
        path: pathlib.Path,
        xyz_data: XyzData,
) -> None:
    xyz_lines = []
    xyz_lines.append(f'{len(xyz_data.elements)}')
    xyz_lines.append('')

    for element, coordinate in zip(xyz_data.elements, xyz_data.coordinates):
        xyz_lines.append(
                f'{element}   {coordinate[0]}   {coordinate[1]}   {coordinate[2]}'
        )
    content = '\n'.join(xyz_lines)
    with open(path, 'w') as xyz_file:
        xyz_file.write(f'{content}\n')


def get_full_xyz_array(
        xyz_data: XyzData,
        removeHs: bool = False,
) -> npt.NDArray:
    if removeHs:
        xyz_array = np.concatenate(
                (np.asarray(xyz_data.elements).T.reshape(-1, 1),
                    xyz_data.coordinates),
                axis=1,
                )
        return np.delete(
                xyz_array,
                np.where(xyz_array[:, 0] == 'H')[0],
                0,
                )
    else:
        return np.concatenate(
            (np.asarray(xyz_data.elements).T.reshape(-1, 1),
                xyz_data.coordinates),
            axis=1,
            )


def get_streusel_radii(
        xyz_data: XyzData,
) -> npt.NDArray[str]:
    streusel_radii = np.array(
            [constants.streusel_radii[element] for element in xyz_data.elements]
    )
    return streusel_radii
