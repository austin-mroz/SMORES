import contextlib
import itertools
import os
import pathlib

import psi4

from smores._internal.esp_molecule import EspMolecule
from smores._internal.molecule import Molecule


def calculate_electrostatic_potential(
    molecule: Molecule | EspMolecule,
    output_directory: pathlib.Path | str,
    grid_dimensions: tuple[int, int, int] = (51, 51, 51),
    num_threads: int | None = None,
    optimize: bool = False,
) -> None:
    """
    Calculate the electrostatic potential.

    Parameters:

        molecule:
            The molecule to optimize.

        outupt_directory:
            The directory in which the calculations are run.

        grid_dimensions:
            The number of voxels in each dimension.

        num_threads:
            The number of threads to use in parallel calculations.
            If ``None``, :func:`os.cpu_count` is used.

        optimize:
            Toggles the optimization of the molecular structure
            before the electrostatic potential is calculated.

    Examples:

        See :mod:`smores.psi4`.

    """

    if num_threads is None:
        num_threads = os.cpu_count()

    output_directory = pathlib.Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)

    _generate_voxel_grid(grid_dimensions, output_directory)

    with _current_working_directory(output_directory):
        psi4.set_options(
            {
                "basis": "aug-cc-pVDZ",
                "CUBEPROP_TASKS": ["ESP"],
                "CUBEPROP_FILEPATH": str(output_directory),
                "reference": "uhf",
            }
        )
        psi4.core.set_num_threads(num_threads)

        psi4_mol = psi4.core.Molecule.from_arrays(
            molecule.positions,
            elem=molecule.atoms,
        )
        psi4.core.set_output_file(
            str(output_directory.joinpath("output.dat")),
            False,
        )

        if optimize:
            psi4.optimize("PBE", molecule=psi4_mol)

        _, wfn = psi4.prop(  # type: ignore
            "PBE",
            molecule=psi4_mol,
            properties=["GRID_ESP"],
            return_wfn=True,
        )
        psi4.cubeprop(wfn)


def _generate_voxel_grid(
    grid_dimensions: tuple[int, int, int],
    output_directory: pathlib.Path,
) -> None:

    x, y, z = grid_dimensions
    # we want to make a box with one corner at (-5,5,5)A and one at
    # (5,5,5)A with a resolution of 0.2 A
    # this should be left to the user eventually
    grid_xyz_coords = []
    for i, j, k in itertools.product(range(x), range(y), range(z)):
        itrans = -5 + 0.2 * i
        jtrans = -5 + 0.2 * j
        ktrans = -5 + 0.2 * k
        grid_xyz_coords.append([itrans, jtrans, ktrans])
    # write the grid to a .dat file
    with open(output_directory.joinpath("grid.dat"), "w") as file:
        for xyz in grid_xyz_coords:
            for c in xyz:
                file.write(str(c) + " ")
            file.write("\n")


@contextlib.contextmanager
def _current_working_directory(directory: pathlib.Path):
    original_directory = os.getcwd()  # aka OGD
    try:
        os.chdir(directory)
        yield
    finally:
        os.chdir(original_directory)
