import pathlib
import typing

import dbstep.Dbstep as db
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

    # def __init__(
    #     self,
    #     atoms: typing.Iterable[str],
    #     positions: npt.ArrayLike,
    #     electrostatic_potential: ElectrostaticPotentialGrid,
    # ) -> None:
    #     """
    #     Initialize an :class:`.EspMolecule`.

    #     Parameters:

    #         atoms (list[str]):
    #             The elemental symbol of each atom of the molecule.

    #         positions (list[list[float]]):
    #             The coordinates of each atom of the molecule,
    #             provided as an N x 3 matrix.

    #         electrostatic_potential:
    #             A 3-D voxel grid of the electrostatic potential
    #             used for calculating the steric parameters.

    #     """

    #     self.atoms = tuple(atoms)
    #     self.positions = np.array(positions)
    #     self._electrostatic_potential = electrostatic_potential

    def get_steric_parameters(
        self,
        dummy_index: int,
        attached_index: int,
    ) -> StericParameters:
        """
        Get the steric parameters from the electrostatic potential.

        Parameters:

            dummy_index:
                The index of the dummy atom.

            attached_index:
                The index of the attached atom of the substituent.

        Returns:
            The parameters.
        """

        db.dbstep(
            atom1=dummy_index,
            atom2=attached_index,
        )

    @classmethod
    def from_cube_file(cls, path: pathlib.Path | str) -> "EspMolecule":
        """
        Get a molecule from a ``.cube`` file.

        Parameters:
            path: The path to the file.
        Returns:
            The molecule.
        """

        instance = cls.__new__(cls)
        instance._cube_file = pathlib.Path(path).resolve()
        return instance
