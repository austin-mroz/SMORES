"""
EspMolecule
===========

"""

import pathlib
import typing

import numpy as np
import numpy.typing as npt

from smores._internal.esp_grid import ElectrostaticPotentialGrid
from smores._internal.steric_parameters import StericParameters


class EspMolecule:
    """
    Calculates :class:`.StericParameters` from electrostatic potentials.

    See Also:

        * :class:`.Molecule`: For calculating steric parameters from
          STREUSEL__ radii.

    __ https://streusel.readthedocs.io

    Examples:

        .. testcode:: get-steric-parameters

            import smores

            molecule = smores.EspMolecule(
                atoms=["H", "Br"],
                positions=[[0., 0., 0.], [1.47, 0., 0.]]
                electrostatic_potential=smores.ElectrostaticPotentialGrid(
                    grid=
                    voxel_size=
                ),
            )
            params = molecule.get_steric_parameters()


    """

    #: The atoms of the molecule.
    atoms: tuple[str, ...]
    #: The N x 3 position matrix of the molecule.
    positions: npt.NDArray[np.float32]

    def __init__(
        self,
        atoms: typing.Iterable[str],
        positions: npt.NDArray[np.float32],
        electrostatic_potential: ElectrostaticPotentialGrid,
    ) -> None:
        """
        Initialize an :class:`.EspMolecule`.

        Parameters:

            atoms:
                The elemental symbol of each atom of the molecule.

            positions:
                The coordinates of each atom of the molecule.

            electrostatic_potential:
                A 3-D voxel grid of the electrostatic potential
                used for calculating the steric parameters.

        """

        self.atoms = tuple(atoms)
        self.positions = np.array(positions)
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
        obj.atoms = cube_data.atoms
        obj.positions = np.array(cube_data.positions)
        obj._electrostatic_potential = cube_data.grid
        return obj

    @classmethod
    def from_mol_file(
        cls,
        path: pathlib.Path | str,
        electrostatic_potential: ElectrostaticPotentialGrid,
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
        obj.atoms = molecule_data.atoms
        obj.positions = np.array(molecule_data.positions)
        obj._electrostatic_potential = electrostatic_potential
        return obj

    @classmethod
    def from_xyz_file(
        cls,
        path: pathlib.Path | str,
        electrostatic_potential: ElectrostaticPotentialGrid,
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
