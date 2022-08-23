import argparse

import numpy as np
from scipy.spatial.distance import cdist
from streusel import gaussian_cube


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "cube_file",
        nargs="+",
        help="The .cube file used to calculate van der Waals volumes.",
    )
    return parser.parse_args()


def _get_distance_matrix(molecule: gaussian_cube.Molecule) -> np.ndarray:
    position_matrix = molecule.coords_and_atoms[["x", "y", "z"]].to_numpy()
    return cdist(position_matrix, position_matrix)


def _get_overlap_matrix(
    distance_matrix: np.ndarray,
    radii_i: np.ndarray,
    radii_j: np.ndarray,
) -> np.ndarray:

    radii_sum = radii_i + radii_j
    overlap = (
        np.pi
        * np.square(radii_sum - distance_matrix)
        * (
            np.square(distance_matrix)
            + 2 * distance_matrix * radii_j
            - 3 * np.square(radii_j)
            + 2 * distance_matrix * radii_i
            - 3 * np.square(radii_i)
            + 6 * radii_j * radii_j
        )
    ) / (12 * distance_matrix)
    overlap[np.isinf(overlap)] = 0
    overlap[radii_sum >= distance_matrix] = 0
    return overlap


def main() -> None:

    vdw_radii = {
        "H": 1.2,
        "C": 1.70,
        "N": 1.55,
        "S": 1.8,
        "O": 1.52,
        "Br": 1.85,
        "Kr": 2.02,
        "Cl": 1.75,
        "Ne": 1.54,
        "F": 1.47,
    }

    cli_args = _get_command_line_arguments()

    for cube_file in cli_args.cube_file:
        mol = gaussian_cube.Molecule(cube_file)
        radii = np.array([vdw_radii[element.symbol] for element in mol.atoms])
        radii_j = np.tile(radii, (len(radii), 1))
        radii_i = radii_j.T
        distance_matrix = _get_distance_matrix(mol)
        overlap_matrix = _get_overlap_matrix(distance_matrix, radii_i, radii_j)
        total_volume = np.sum(4 * np.pi * np.power(radii, 3) / 3)
        total_overlap = np.sum(overlap_matrix) / 2
        print(
            f"{cube_file} "
            f"-- VOLUME: {total_volume} "
            f"-- OVERLAP: {total_overlap} "
            f"-- DIFFERENCE: {total_volume - total_overlap}"
        )


if __name__ == "__main__":
    main()
