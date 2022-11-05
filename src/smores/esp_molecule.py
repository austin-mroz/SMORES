"""
EspMolecule
===========

"""

import pathlib
import typing
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

from smores._internal import file_readers
from smores.steric_parameters import StericParameters


@dataclass(frozen=True, slots=True)
class ElectrosaticPotentialGrid:
    """
    A 3-D grid of voxels holding electrostatic potentials.

    Attributes:
        grid: The voxel grid.
        voxl_size: The length of a single voxel in each dimension.

    """

    grid: npt.ArrayLike
    voxel_size: tuple[float, float, float]


class EspMolecule:
    """
    Calculates :class:`.StericParameters` from electrostatic potentials.

    Examples:

        .. testcode:: get-steric-parameters

            import smores

            molecule = smores.EspMolecule(
                atoms=["H", "Br"],
                positions=
                electrostatic_potential=smores.ElectrosaticPotentialGrid(
                    grid=
                    voxel_size=
                ),
            )
            molecule.get_steric_parameters()

    See Also:

        * :class:`.Molecule`

    """

    def __init__(
        self,
        atoms: typing.Iterable[str],
        positions: npt.ArrayLike,
        electrostatic_potential: ElectrosaticPotentialGrid,
    ) -> None:
        """
        Initialize a :class:`.EspMolecule`.

        Parameters:

            atoms:
                The elemental symbol of each atom of the molecule.

            positions:
                The coordinates of each atom of the molecule.

            electrostatic_potential:
                A 3-D voxel grid of the electrostatic potential
                used for calculating the steric parameters.

        """

        self._atoms = tuple(atoms)
        self._postions = np.array(positions)
        self._electrostatic_potential = electrostatic_potential

    def get_steric_parameters(self) -> StericParameters:
        """
        Get the steric parameters from the electrostatic potential.

        Returns:
            The parameters.
        """

        pass

    @classmethod
    def from_cube_file(cls, path: pathlib.Path | str) -> "EspMolecule":
        """
        Get a molecule from a ``.cube`` file.

        Parameters:
            path: The path to the file.
        Returns:
            The molecule.
        """

        path = pathlib.Path(path)
        cube_data = file_readers.read_cube_file(path)
        obj = cls.__new__(cls)
        obj._atoms = cube_data.atoms
        obj._postions = np.array(cube_data.positions)
        obj._electrostatic_potential = cube_data.grid
        return obj

    @classmethod
    def from_mol_file(
        cls,
        path: pathlib.Path | str,
        electrostatic_potential: ElectrosaticPotentialGrid,
    ) -> "EspMolecule":
        """
        Get a molecule from a ``.mol`` file.

        Parameters:
            path: The path to the file.
        Returns:
            The molecule.
        """

        path = pathlib.Path(path)
        molecule_data = file_readers.read_mol_file(path)
        obj = cls.__new__(cls)
        obj._atoms = molecule_data.atoms
        obj._postions = np.array(molecule_data.positions)
        obj._electrostatic_potential = electrostatic_potential
        return obj

    @classmethod
    def from_xyz_file(
        cls,
        path: pathlib.Path | str,
        electrostatic_potential: ElectrosaticPotentialGrid,
    ) -> "EspMolecule":
        """
        Get a molecule from a ``.xyz`` file.

        Parameters:
            path: The path to the file.
        Returns:
            The molecule.
        """

        path = pathlib.Path(path)
        return cls()
