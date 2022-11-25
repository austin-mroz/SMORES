import pathlib
import typing

import morfeus
import dbstep.Dbstep as db
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

    @classmethod
    def from_xyz_file(
        cls,
        path: pathlib.Path | str,
        radii: npt.ArrayLike | None = None,
    ) -> "Molecule":
        """
        Get a molecule from a ``.xyz`` file.

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
        conformer_id=0,
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

    def get_steric_parameters_using_radii(
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

    def get_steric_parameters_using_cube_file(
            self,
            dummy_index: int,
            attached_index: int,
            cube_file: str | pathlib.Path,
    ) -> StericParameters:
        """
        Get the steric parameters from STREUSEL_ cube file.
        
        Parameters:
        
            dummy_index:
                The index of the dummy atom.

            attached_index:
                The index of the attached atom of the substituent.

        Returns:
            The steric parameters.

        .. _STREUSEL: https://streusel.readthedocs.io
        
        """

        sterimol = db.dbstep(
                str(cube_file),
                atom1=dummy_index,
                atom2=attached_index,
                sterimol=True,
                surface="density",
        )

        return StericParameters(
                L=sterimol.L,
                B1=sterimol.Bmin,
                B5=sterimol.Bmax,
        )
