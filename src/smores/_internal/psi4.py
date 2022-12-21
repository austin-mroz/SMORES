import itertools
import os
import pathlib
import typing

import numpy as np
import numpy.typing as npt
import psi4
import rdkit.Chem.AllChem as rdkit

from smores._internal.constants import atomic_mass
from smores._internal.utilities import current_working_directory


def calculate_electrostatic_potential(
    molecule: rdkit.Mol,
    output_directory: pathlib.Path | str,
    grid_origin: tuple[float, float, float],
    grid_length: float,
    num_voxels_per_dimension: int,
    num_threads: int | None = None,
    optimize: bool = False,
    conformer_id: int = 0,
) -> pathlib.Path:
    """
    Calculate the electrostatic potential.

    Parameters:

        molecule:
            The molecule to optimize. Must have at
            least one conformer.

        output_directory:
            The directory in which the calculations are run.

        grid_origin:
            The origin of the grid.

        grid_length:
            The length of the grid in each dimension in Angstrom.

        num_voxels_per_dimension:
            The number of voxels in each dimension.

        num_threads:
            The number of threads to use in parallel calculations.
            If ``None``, :func:`os.cpu_count` is used.

        optimize:
            Toggles the optimization of the molecular structure
            before the electrostatic potential is calculated.

        conformer_id:
            The conformer of `molecule` to use.

    Returns:
        The path to the ``.cube`` file holding the electrostatic
        potential of `molecule`.

    Examples:

        *Calculate the steric parameters from a cube file*

        .. testcode:: calculate-steric-parameters-from-cube-file

            import smores
            import smores.psi4

            molecule = smores.rdkit_from_smiles("Br")
            cube_path = smores.psi4.calculate_electrostatic_potential(
                molecule=molecule,
                output_directory="HBr_ESP",
                grid_origin=(-10, -10, -10),
                grid_length=20,
                num_voxels_per_dimension=50,
            )
            esp_molecule = smores.EspMolecule.from_cube_file(cube_path, \
dummy_index=0, attached_index=1)

        .. doctest:: calculate-steric-parameters-from-cube-file

            >>> esp_molecule.get_steric_parameters()
            StericParameters(L=3.57164113574581, \
B1=1.9730970556668774, B5=2.320611610648539)

    """

    if num_threads is None:
        num_threads = os.cpu_count()

    output_directory = pathlib.Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)

    _generate_voxel_grid(
        output_directory=output_directory,
        grid_origin=grid_origin,
        grid_length=grid_length,
        num_voxels_per_dimension=num_voxels_per_dimension,
    )

    with current_working_directory(output_directory):
        psi4.set_options(
            {
                "basis": "aug-cc-pVDZ",
                "CUBEPROP_TASKS": ["ESP"],
                "CUBEPROP_FILEPATH": str(output_directory),
                "reference": "uhf",
            }
        )
        psi4.core.set_num_threads(num_threads)

        elements = [atom.GetSymbol() for atom in molecule.GetAtoms()]
        coordinates = molecule.GetConformer(conformer_id).GetPositions()
        center_of_mass = _get_center_of_mass(
            elements=elements,
            coordinates=coordinates,
        )
        psi4_mol = psi4.core.Molecule.from_arrays(
            coordinates - center_of_mass,
            elem=elements,
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
        return output_directory / "ESP.cube"


def _generate_voxel_grid(
    output_directory: pathlib.Path,
    grid_origin: tuple[float, float, float],
    grid_length: float,
    num_voxels_per_dimension: int,
) -> None:

    origin_x, origin_y, origin_z = grid_origin
    voxel_size = grid_length / num_voxels_per_dimension

    grid_xyz_coords = []
    for i, j, k in itertools.product(
        range(num_voxels_per_dimension),
        range(num_voxels_per_dimension),
        range(num_voxels_per_dimension),
    ):
        itrans = -origin_x + voxel_size * i
        jtrans = -origin_y + voxel_size * j
        ktrans = -origin_z + voxel_size * k
        grid_xyz_coords.append([itrans, jtrans, ktrans])

    with open(output_directory.joinpath("grid.dat"), "w") as file:
        for xyz in grid_xyz_coords:
            for c in xyz:
                file.write(str(c) + " ")
            file.write("\n")


def _get_center_of_mass(
    elements: typing.Iterable[str],
    coordinates: npt.NDArray[np.float32],
) -> npt.NDArray[np.float32]:

    atom_masses = np.array([atomic_mass[element] for element in elements])
    scaled_coordinates = coordinates * atom_masses[:, np.newaxis]
    return scaled_coordinates.sum(axis=0) / atom_masses.sum()
