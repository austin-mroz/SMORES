import pathlib
import tempfile
import typing

import ase
import ase.io.cube
import dbstep.Dbstep as db
import numpy as np
import numpy.typing as npt
import rdkit.Chem.AllChem as rdkit
import streusel.gaussian_cube

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

        """

        self._atoms = tuple(atoms)
        self._positions = np.array(positions)
        self._dummy_index = dummy_index
        self._attached_index = attached_index

        with tempfile.TemporaryDirectory() as tmp_dir:
            electrostatic_potential_file = pathlib.Path(tmp_dir) / "esp.cube"
            write_cube(
                path=electrostatic_potential_file,
                voxels=electrostatic_potential.voxels,
                positions=self._positions,
                elements=self._atoms,
                voxel_origin=electrostatic_potential.voxel_origin,
                voxel_dimensions=electrostatic_potential.voxel_size,
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
            voxel_size=electrostatic_potential.voxel_size,
            voxel_origin=electrostatic_potential.voxel_origin,
        )

    def get_steric_parameters(self) -> StericParameters:
        """
        Get the steric parameters from the electrostatic potential.

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
                voxel_dimensions=self._electric_field_surface.voxel_size,
            )

            params = db.dbstep(
                str(electric_field_surface_file),
                atom1=self._dummy_index,
                atom2=self._attached_index,
                surface="density",
                sterimol=True,
                quiet=True,
            )

        return StericParameters(
            L=params.L,
            B1=params.Bmin,
            B5=params.Bmax,
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
        grid, atoms = ase.io.cube.read_cube_data(str(path))
        instance._atoms = atoms.get_atomic_numbers()
        instance._positions = atoms.positions / ase.units.Bohr
        instance._dummy_index = dummy_index
        instance._attached_index = attached_index

        streusel_molecule = streusel.gaussian_cube.Molecule(str(path))

        instance._electric_field_surface = VoxelGrid(
            voxels=_get_electric_field_surface(streusel_molecule),
            voxel_size=streusel_molecule.res,
            voxel_origin=streusel_molecule.origin,
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


def _get_electric_field_surface(
    streusel_molecule: streusel.gaussian_cube.Molecule,
) -> npt.NDArray:
    streusel_molecule.get_efield()
    streusel_molecule.sample_efield()
    return streusel_molecule.surface_mask
