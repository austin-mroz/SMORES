import math
import pathlib
import tempfile
import typing

import dbstep.Dbstep as db
import flour
import numpy as np
import numpy.typing as npt
import rdkit.Chem.AllChem as rdkit
import streusel.gaussian_cube

from smores._internal.from_esp import calculate_steric_parameters_from_esp
from smores._internal.steric_parameters import StericParameters
from smores._internal.voxel_grid import VoxelGrid


class EspMolecule:
    """
    Calculates :class:`.StericParameters` from electrostatic potentials.

    See Also:

        * :class:`.Molecule`: For calculating steric parameters from
          STREUSEL__ radii.

    __ https://streusel.readthedocs.io

    Examples:

        *Calculate steric parameters*

        .. doctest:: esp-molecule-calculate-steric-parameters

            >>> import smores
            >>> molecule = smores.EspMolecule.from_cube_file("HBr", \
dummy_index=0, attached_index=1)
            >>> molecule.get_steric_parameters()
            StericParameters(L=3.57164113574581, \
B1=1.9730970556668774, B5=2.320611610648539)


    """

    def __init__(
        self,
        atoms: typing.Iterable[str],
        positions: npt.ArrayLike,
        dummy_index: int,
        attached_index: int,
        electrostatic_potential: VoxelGrid,
        excluded_indices: typing.Iterable[int] | None = None,
    ) -> None:
        """
        Initialize an :class:`.EspMolecule`.

        Parameters:

            atoms (list[str]):
                The elemental symbol of each atom of the molecule.

            positions (list[list[float]]):
                The coordinates of each atom of the molecule,
                provided as an N x 3 matrix.

            dummy_index:
                The index of the dummy atom.

            attached_index:
                The index of the attached atom of the substituent.

            electrostatic_potential:
                The electrostatic potential used for calculating
                the steric parameters.

            excluded_indices (list[int]):
                The indices of atoms which are not included in the
                parameter calculation.

        """

        self._atoms = np.array(
            [rdkit.Atom(atom).GetAtomicNum() for atom in atoms],
            dtype=np.uint8,
        )
        self._positions = np.array(positions, dtype=np.float64)
        self._dummy_index = dummy_index
        self._attached_index = attached_index

        with tempfile.TemporaryDirectory() as tmp_dir:
            electrostatic_potential_file = pathlib.Path(tmp_dir) / "esp.cube"
            flour.write_cube(
                path=electrostatic_potential_file,
                title1="t1",
                title2="t2",
                atoms=self._atoms,
                charges=np.zeros(len(self._atoms), dtype=np.float64),
                positions=self._positions,
                voxel_origin=electrostatic_potential.voxel_origin,
                voxel_size=electrostatic_potential.voxel_size,
                voxels=electrostatic_potential.voxels,
            )
            streusel_molecule = streusel.gaussian_cube.Molecule(
                electrostatic_potential_file,
            )
            streusel_molecule.get_efield()
            streusel_molecule.sample_efield()
            electric_field_surface = np.zeros(
                electrostatic_potential.voxels.shape
            )
            electric_field_surface[streusel_molecule.surface_ijk] = 1

        if excluded_indices is not None:
            exclude = np.ones(electrostatic_potential.voxels.shape)
            exclude_size = 10
            for index in excluded_indices:
                position = self._positions[index]
                voxel = position // electrostatic_potential.voxel_size
                voxel_range = slice(voxel - exclude_size, voxel + exclude_size)
                exclude[voxel_range, voxel_range, voxel_range] = False
            electric_field_surface[exclude] = False
        self._electric_field_surface = VoxelGrid(
            voxels=electric_field_surface,
            voxel_size=electrostatic_potential.voxel_size,
            voxel_origin=electrostatic_potential.voxel_origin,
        )

    def get_electric_field_surface(self) -> VoxelGrid:
        return VoxelGrid(
            voxels=np.array(self._electric_field_surface.voxels),
            voxel_size=np.array(self._electric_field_surface.voxel_size),
            voxel_origin=np.array(self._electric_field_surface.voxel_origin),
        )

    def get_steric_parameters(self) -> StericParameters:
        """
        Get the steric parameters from the electrostatic potential.

        Returns:
            The parameters.

        """
        cube_data_positions_idx = _convert_euclidean_positions_to_indices(
            positions=self._positions,
            origin=self._electric_field_surface.voxel_origin,
            resolution=self._electric_field_surface.voxel_size,
            voxel_grid_shape=self._electric_field_surface.voxels.shape,
        )
        return calculate_steric_parameters_from_esp(
            self._electric_field_surface.voxels.astype(float),
            self._electric_field_surface.voxel_size * np.identity(3),
            cube_data_positions_idx[self._attached_index],
            cube_data_positions_idx[self._dummy_index],
        )

    def get_dummy_index(self) -> int:
        """
        Get the index of the dummy atom.

        Returns:
            The index of the dummy atom.

        """
        return self._dummy_index

    def get_attached_index(self) -> int:
        """
        Get the index of the atom attached to the substituent.

        Returns:
            The index of the atom attached to the substituent.

        """
        return self._attached_index

    @classmethod
    def from_cube_file(
        cls,
        path: pathlib.Path | str,
        dummy_index: int,
        attached_index: int,
    ) -> "EspMolecule":
        """
        Get a molecule from a ``.cube`` file.

        Parameters:

            path:
                The path to the file.

            dummy_index:
                The index of the dummy atom.

            attached_index:
                The index of the attached atom of the substituent.

        Returns:
            The molecule.

        """

        instance = cls.__new__(cls)
        cube_data = flour.read_cube(path)
        instance._atoms = cube_data.atoms
        instance._positions = cube_data.positions
        instance._dummy_index = dummy_index
        instance._attached_index = attached_index

        streusel_molecule = streusel.gaussian_cube.Molecule(str(path))

        instance._electric_field_surface = VoxelGrid(
            voxels=_get_electric_field_surface(streusel_molecule),
            voxel_size=cube_data.grid.voxel_size.sum(axis=0),
            voxel_origin=cube_data.grid.origin,
        )

        return instance

    @classmethod
    def from_rdkit(
        cls,
        molecule: rdkit.Mol,
        dummy_index: int,
        attached_index: int,
        electrostatic_potential: VoxelGrid,
        conformer_id: int = 0,
    ) -> "EspMolecule":
        """
        Get a molecule from an :mod:`rdkit` molecule.

        Parameters:

            molecule:
                The :mod:`rdkit` molecule. It must have at least
                one conformer.

            dummy_index:
                The index of the dummy atom.

            attached_index:
                The index of the attached atom of the substituent.

            electrostatic_potential:
                The electrostatic potential used for calculating
                the steric parameters.

            conformer_id:
                The id of the conformer to use.

        Returns:
            The :mod:`smores` molecule.

        """

        return cls(
            atoms=(atom.GetSymbol() for atom in molecule.GetAtoms()),
            positions=molecule.GetConformer(conformer_id).GetPositions(),
            dummy_index=dummy_index,
            attached_index=attached_index,
            electrostatic_potential=electrostatic_potential,
        )


def _convert_euclidean_positions_to_indices(
    positions,
    origin,
    resolution,
    voxel_grid_shape,  # : flour.CubeData,
) -> npt.NDArray:
    # translate origin and positions to 0,0,0
    translated_position_matrix = positions - origin

    lattice_lengths = resolution.sum(axis=0) * np.asarray(voxel_grid_shape)

    percent_along_lattice_vector = np.divide(
        translated_position_matrix, lattice_lengths
    )

    return (percent_along_lattice_vector * voxel_grid_shape).astype(int)


def _get_electric_field_surface(
    streusel_molecule: streusel.gaussian_cube.Molecule,
) -> npt.NDArray:
    streusel_molecule.get_efield()
    streusel_molecule.sample_efield()
    return streusel_molecule.surface_mask
