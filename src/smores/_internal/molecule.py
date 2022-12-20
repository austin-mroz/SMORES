import pathlib
import typing

import morfeus
import numpy as np
import numpy.typing as npt
import rdkit.Chem.AllChem as rdkit

from smores._internal.constants import streusel_radii
from smores._internal.steric_parameters import StericParameters


class Molecule:
    """
    Calculates :class:`.StericParameters` from STREUSEL radii.

    See also:

        * :class:`.EspMolecule`: For calculating steric parameters
          from electrostatic potentials.

    Examples:

        Using custom atomic radii

        .. testcode:: custom-atomic-radii

            import smores
            molecule = smores.Molecule.from_smiles("CBr")
            params = molecule.get_steric_parameters(0, 1)

    """

    #: The atoms of the molecule.
    _atoms: tuple[str, ...]
    #: The N x 3 position matrix of the molecule.
    _positions: npt.NDArray[np.float32]

    def __init__(
        self,
        atoms: typing.Iterable[str],
        positions: npt.ArrayLike,
        dummy_index: int,
        attached_index: int,
        excluded_indices: typing.Iterable[int] | None = None,
        radii: npt.ArrayLike | None = None,
    ) -> None:
        """
        Initialize a :class:`.Molecule`.

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

            excluded_indices (list[int]):
                The indices of atoms which are not included in the
                parameter calculation.

            radii (list[float]):
                The radius of each atom of the molecule. If
                ``None`` the STREUSEL_ radii will be used.

        .. _STREUSEL: https://streusel.readthedocs.io


        """

        if radii is None:
            radii = np.array([streusel_radii[atom] for atom in atoms])
        else:
            radii = np.array(radii)

        self._atoms = tuple(atoms)
        self._positions = np.array(positions)
        self._radii = radii
        self._morfeus_dummy_index = dummy_index + 1
        self._morfeus_attached_index = attached_index + 1
        self._morfeus_excluded_indices = (
            None
            if excluded_indices is None
            else tuple(index + 1 for index in excluded_indices)
        )

    def get_dummy_index(self) -> int:
        return self._morfeus_dummy_index - 1

    def get_attached_index(self) -> int:
        return self._morfeus_attached_index - 1

    def get_excluded_indices(self) -> typing.Iterable[int]:
        if self._morfeus_excluded_indices is None:
            return
        for index in self._morfeus_excluded_indices:
            yield index - 1

    @classmethod
    def from_xyz_file(
        cls,
        path: pathlib.Path | str,
        dummy_index: int,
        attached_index: int,
        excluded_indices: typing.Iterable[int] | None = None,
        radii: npt.ArrayLike | None = None,
    ) -> "Molecule":
        """
        Get a molecule from a ``.xyz`` file.

        Parameters:

            path:
                The path to the file.

            dummy_index:
                The index of the dummy atom.

            attached_index:
                The index of the attached atom of the substituent.

            excluded_indices (list[int]):
                The indices of atoms which are not included in the
                parameter calculation.

            radii (list[float]):
                The radius of each atom of the molecule. If
                ``None`` the STREUSEL_ radii will be used.

        Returns:
            The molecule.

        """

        instance = cls.__new__(cls)
        molecule = rdkit.MolFromXYZFile(str(path))
        instance._atoms = tuple(
            atom.GetSymbol() for atom in molecule.GetAtoms()
        )
        instance._positions = molecule.GetConformer(0).GetPositions()
        if radii is None:
            instance._radii = np.array(
                [streusel_radii[atom] for atom in instance._atoms]
            )
        else:
            instance._radii = np.array(radii)
        instance._morfeus_dummy_index = dummy_index + 1
        instance._morfeus_attached_index = attached_index + 1
        instance._morfeus_excluded_indices = (
            None
            if excluded_indices is None
            else tuple(index + 1 for index in excluded_indices)
        )
        return instance

    @classmethod
    def from_mol_file(
        cls,
        path: pathlib.Path | str,
        radii: npt.ArrayLike | None = None,
    ) -> "Molecule":
        """
        Get a molecule from a ``.mol`` file.

        Parameters:

            path:
                The path to the file.

            radii (list[float]):
                The radius of each atom of the molecule. If
                ``None`` the STREUSEL_ radii will be used.

        Returns:
            The molecule.
        """

        instance = cls.__new__(cls)
        molecule = rdkit.MolFromMolFile(str(path), removeHs=False)
        instance._atoms = tuple(
            atom.GetSymbol() for atom in molecule.GetAtoms()
        )
        instance._positions = molecule.GetConformer(0).GetPositions()
        if radii is None:
            instance._radii = np.array(
                [streusel_radii[atom] for atom in instance._atoms]
            )
        else:
            instance._radii = np.array(radii)
        return instance

    @classmethod
    def from_smiles(
        cls,
        smiles: str,
        positions: npt.ArrayLike | None = None,
        radii: npt.ArrayLike | None = None,
    ) -> "Molecule":
        """
        Get a molecule from a SMILES string.

        Parameters:

            smiles:
                The SMILES of the molecule.

            positions (list[list[float]]):
                The coordinates of each atom of the molecule
                provided as an N x 3 matrix.
                If ``None`` then the molecule will have its
                coordinates calculated with ETKDG_.

            radii (list[float]):
                The radius of each atom of the molecule. If
                ``None`` the STREUSEL_ radii will be used.

        Returns:
            The molecule.

        .. _ETKDG: https://www.rdkit.org/docs/source/\
rdkit.Chem.rdDistGeom.html#rdkit.Chem.rdDistGeom.ETKDGv3

        """

        instance = cls.__new__(cls)
        molecule = rdkit.AddHs(rdkit.MolFromSmiles(smiles))
        instance._atoms = tuple(
            atom.GetSymbol() for atom in molecule.GetAtoms()
        )

        if positions is None:
            params = rdkit.ETKDGv3()
            params.randomSeed = 4
            rdkit.EmbedMolecule(molecule, params)
            instance._positions = molecule.GetConformer(0).GetPositions()
        else:
            instance._positions = np.array(positions)

        if radii is None:
            instance._radii = np.array(
                [streusel_radii[atom] for atom in instance._atoms]
            )
        else:
            instance._radii = np.array(radii)

        return instance

    @classmethod
    def from_rdkit(
        cls,
        molecule: rdkit.Mol,
        radii: npt.ArrayLike | None = None,
        conformer_id: int = 0,
    ) -> "Molecule":
        """
        Get a molecule from an :mod:`rdkit` molecule.

        Parameters:

            molecule:
                The :mod:`rdkit` molecule. It must have at least
                one conformer.

            radii (list[float]):
                The radius of each atom of the molecule. If
                ``None`` the STREUSEL_ radii will be used.

            conformer_id:
                The id of the conformer to use.

        Returns:
            The :mod:`smores` molecule.


        """

        instance = cls.__new__(cls)
        instance._atoms = tuple(
            atom.GetSymbol() for atom in molecule.GetAtoms()
        )
        instance._positions = molecule.GetConformer(
            conformer_id
        ).GetPositions()

        if radii is None:
            instance._radii = np.array(
                [streusel_radii[atom] for atom in instance._atoms]
            )
        else:
            instance._radii = np.array(radii)

        return instance

    def get_steric_parameters(self) -> StericParameters:
        """
        Get the steric paramters from STREUSEL_ radii.

        Returns:
            The parameters.

        .. _STREUSEL: https://streusel.readthedocs.io

        """

        sterimol = morfeus.Sterimol(
            elements=self._atoms,
            coordinates=self._positions,
            dummy_index=self._morfeus_dummy_index,
            attached_index=self._morfeus_attached_index,
            radii=self._radii,
            excluded_atoms=self._morfeus_excluded_indices,
        )
        return StericParameters(
            L=sterimol.L_value,
            B1=sterimol.B_1_value,
            B5=sterimol.B_5_value,
        )
