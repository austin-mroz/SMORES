import os
import pathlib
import typing
from collections import abc
from itertools import product

import morfeus
import numpy as np
import numpy.typing as npt
import psi4
import rdkit.Chem.AllChem as rdkit

from smores._internal.constants import streusel_radii
from smores._internal.steric_parameters import StericParameters


class InvalidDirectoryError(Exception):
    pass


class Molecule:
    """
    Calculates :class:`.StericParameters` from `STREUSEL`_ radii.


    .. _STREUSEL: https://streusel.readthedocs.io

    See also:

        * :class:`.EspMolecule`: For calculating steric parameters
          from electrostatic potentials.

    Examples:

        Using custom atomic radii

        .. testcode:: custom-atomic-radii

            import smores
            molecule = smores.Molecule.from_smiles("CBr")
            params = molecule.get_steric_parameters(
                dummy_index=0,
                attached_index=1,
            )



    """

    _atoms: abc.Sequence[str]

    def __init__(
        self,
        atoms: typing.Iterable[str],
        positions: npt.NDArray[np.float32],
        radii: npt.NDArray[np.float32] | None = None,
    ) -> None:
        """
        Initialize a :class:`.Molecule`.

        Parameters:

            atoms:
                The elemental symbol of each atom of the molecule.

            positions:
                The coordinates of each atom of the molecule.

            radii:
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

    @classmethod
    def from_xyz_file(
        cls,
        path: pathlib.Path | str,
        radii: npt.NDArray[np.float32] | None = None,
    ) -> "Molecule":
        """
        Get a molecule from a ``.xyz`` file.

        Parameters:

            path:
                The path to the file.

            radii:
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
        return instance

    @classmethod
    def from_mol_file(
        cls,
        path: pathlib.Path | str,
        radii: npt.NDArray[np.float32] | None = None,
    ) -> "Molecule":
        """
        Get a molecule from a ``.mol`` file.

        Parameters:

            path:
                The path to the file.

            radii:
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
        positions: npt.NDArray[np.float32] | None = None,
        radii: npt.NDArray[np.float32] | None = None,
    ) -> "Molecule":
        """
        Get a molecule from a SMILES string.

        Parameters:

            smiles:
                The SMILES of the molecule.

            positions:
                The coordinates of each atom of the moleclue.
                If ``None`` then the molecule will have its
                coordinates calculated with ETKDG_.

            radii:
                The radius of each atom of the molecule. If
                ``None`` the STREUSEL_ radii will be used.

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

    def get_steric_parameters(
        self,
        dummy_index: int,
        attached_index: int,
    ) -> StericParameters:
        """
        Get the steric paramters from STREUSEL_ radii.


        Parameters:

            dummy_index:
                The index of the dummy atom.

            attached_index:
                The index of the attached atom of the substituent.

        Returns:
            The parameters.

        .. _STREUSEL: https://streusel.readthedocs.io

        """

        sterimol = morfeus.Sterimol(
            elements=self._atoms,
            coordinates=self._positions,
            dummy_index=dummy_index + 1,
            attached_index=attached_index + 1,
            radii=self._radii,
        )
        return StericParameters(
            L=sterimol.L_value,
            B1=sterimol.B_1_value,
            B5=sterimol.B_5_value,
        )
