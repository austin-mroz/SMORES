import pathlib
import tempfile
import typing
from collections import abc

import dbstep.Dbstep as db
import numpy as np
import numpy.typing as npt
import streusel.gaussian_cube

from smores._internal.read_cube import read_cube
from smores._internal.steric_parameters import StericParameters
from smores._internal.voxel_grid import VoxelGrid
from smores._internal.write_cube import write_cube


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

    def __init__(
        self,
        atoms: typing.Iterable[int],
        positions: npt.ArrayLike,
        electrostatic_potential: VoxelGrid,
    ) -> None:
        """
        Initialize an :class:`.EspMolecule`.

        Parameters:

            atoms (list[int]):
                The atomic number of each atom of the molecule.

            positions (list[list[float]]):
                The coordinates of each atom of the molecule,
                provided as an N x 3 matrix.

            electrostatic_potential:
                The electrostatic potential used for calculating
                the steric parameters.

        """

        self._atoms: abc.Collection = tuple(atoms)
        self._positions = np.array(positions)

        with tempfile.TemporaryDirectory() as tmp_dir:
            electrostatic_potential_file = pathlib.Path(tmp_dir) / "esp.cube"
            write_cube(
                path=electrostatic_potential_file,
                voxels=electrostatic_potential.voxels,
                positions=self._positions,
                elements=self._atoms,
                voxel_origin=electrostatic_potential.voxel_origin,
                voxel_x_vector=electrostatic_potential.voxel_x_vector,
                voxel_y_vector=electrostatic_potential.voxel_y_vector,
                voxel_z_vector=electrostatic_potential.voxel_z_vector,
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

        self._electric_field_surface = VoxelGrid(
            voxels=electric_field_surface,
            voxel_origin=electrostatic_potential.voxel_origin,
            voxel_x_vector=electrostatic_potential.voxel_x_vector,
            voxel_y_vector=electrostatic_potential.voxel_y_vector,
            voxel_z_vector=electrostatic_potential.voxel_z_vector,
        )

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

        with tempfile.TemporaryDirectory() as tmp_dir:
            electric_field_surface_file = (
                pathlib.Path(tmp_dir) / "ef_surface.cube"
            )
            write_cube(
                path=electric_field_surface_file,
                voxels=self._electric_field_surface.voxels,
                positions=self._positions,
                elements=self._atoms,
                voxel_origin=self._electric_field_surface.voxel_origin,
                voxel_x_vector=self._electric_field_surface.voxel_x_vector,
                voxel_y_vector=self._electric_field_surface.voxel_y_vector,
                voxel_z_vector=self._electric_field_surface.voxel_z_vector,
            )

            params = db.dbstep(
                str(electric_field_surface_file),
                atom1=dummy_index,
                atom2=attached_index,
                surface="density",
                sterimol=True,
                quiet=True,
            )

        return StericParameters(
            L=params.L,
            B1=params.Bmin,
            B5=params.Bmax,
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
        cube_data = read_cube(pathlib.Path(path))
        instance._atoms = cube_data.atomic_numbers
        instance._positions = np.array(cube_data.positions)
        streusel_molecule = streusel.gaussian_cube.Molecule(str(path))

        instance._electric_field_surface = VoxelGrid(
            voxels=_get_electric_field_surface(streusel_molecule),
            voxel_origin=cube_data.voxel_grid.voxel_origin,
            voxel_x_vector=cube_data.voxel_grid.voxel_x_vector,
            voxel_y_vector=cube_data.voxel_grid.voxel_y_vector,
            voxel_z_vector=cube_data.voxel_grid.voxel_z_vector,
        )

        return instance


def _get_electric_field_surface(
    streusel_molecule: streusel.gaussian_cube.Molecule,
) -> npt.NDArray:
    streusel_molecule.get_efield()
    streusel_molecule.sample_efield()
    return streusel_molecule.surface_mask
