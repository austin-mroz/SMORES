import numpy as np
import itertools

def write_cube(
    filepath:str,
    vox_grid:np.ndarray,
    coords:np.ndarray,
    elements:list[str],
    origin: np.ndarray,
    res: np.ndarray,
) -> None:
    origin_x, origin_y, origin_z = origin
    res_x, res_y, res_z = res
    nvox_x, nvox_y, nvox_z = vox_grid.shape
    with open(filepath, 'w') as cube:
        # write header stuff
        lines = [
            ' title',
            ' title2',
            f'{len(elements): >5} {origin_x: >11.6f} {origin_y: >11.6f} {origin_z: >11.6f}',
            f'{nvox_x: >5} {res_x: >11.6f} {0.: >11.6f} {0.: >11.6f}',
            f'{nvox_y: >5} {0.: >11.6f} {res_y: >11.6f} {0.: >11.6f}',
            f'{nvox_z: >5} {0.: >11.6f} {0.: >11.6f} {res_z: >11.6f}',
        ]
        for element, [x, y, z] in zip(elements, coords):
            lines.append(f'{element: >5} {0.: >11.6f} {x: >11.6f} {y: >11.6f} {z: >11.6f}')

        # write voxel grid to cube file from tensor
        ncols = 6
        for i in range(nvox_x):
            for j in range(nvox_y):
                line = []
                for k in range(nvox_z):
                    line.append(f'{vox_grid[i, j, k]: >12.5E}')
                    if len(line) == ncols:
                        lines.append(_pad_left(" ".join(line)))
                        line = []
                lines.append(_pad_left(" ".join(line)))
                line = []
        content = '\n'.join(lines)
        cube.write(f'{content}\n')


def _pad_left(s: str) -> str:
    return f" {s}"
