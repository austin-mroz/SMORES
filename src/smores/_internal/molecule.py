import os
import pathlib
import typing
from collections import abc
from itertools import product

import numpy as np
import numpy.typing as npt
import psi4
import rdkit.Chem.AllChem as rdkit

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
            params = molecule.get_steric_parameters({"C": 1.6})



    """

    _atoms: abc.Sequence[str]

    def __init__(
        self,
        atoms: typing.Iterable[str],
        positions: npt.NDArray[np.float32],
    ) -> None:
        """
        Initialize a :class:`.Molecule`.

        Parameters:

            atoms:
                The elemental symbol of each atom of the molecule.

            positions:
                The coordinates of each atom of the molecule.

        """

        self._atoms = tuple(atoms)
        self._positions = np.array(positions)

    @classmethod
    def from_xyz_file(
        cls,
        path: pathlib.Path | str,
    ) -> "Molecule":
        """
        Get a molecule from a ``.xyz`` file.

        Parameters:
            path: The path to the file.
        Returns:
            The molecule.

        """

        instance = cls.__new__(cls)
        molecule = rdkit.MolFromXYZFile(str(path))
        instance._atoms = tuple(
            atom.GetSymbol() for atom in molecule.GetAtoms()
        )
        instance._positions = molecule.GetConformer(0).GetPositions()
        return instance

    @classmethod
    def from_mol_file(
        cls,
        path: pathlib.Path | str,
    ) -> "Molecule":
        """
        Get a molecule from a ``.mol`` file.

        Parameters:
            path: The path to the file.
        Returns:
            The molecule.
        """

        instance = cls.__new__(cls)
        molecule = rdkit.MolFromMolFile(str(path))
        instance._atoms = tuple(
            atom.GetSymbol() for atom in molecule.GetAtoms()
        )
        instance._positions = molecule.GetConformer(0).GetPositions()
        return instance

    @classmethod
    def from_smiles(
        cls,
        smiles: str,
        positions: npt.NDArray[np.float32] | None = None,
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
        return instance

    def get_steric_parameters(
        self,
        radii: abc.Mapping[str, float] | None = None,
    ) -> StericParameters:
        """
        Get the steric paramters from STREUSEL_ radii.


        Parameters:
            radii:
                The atomic radii to use when calculating the steric
                parameters. If ``None`` then standard STREUSEL
                atomic radii will be used.

        Returns:
            The parameters.

        .. _STREUSEL: https://streusel.readthedocs.io

        """

        pass

    def _generate_voxel_grid(
        self,
        resolution: int,
        output_directory: pathlib.Path,
    ) -> None:
        # we want to make a box with one corner at (-5,5,5)A and one at
        # (5,5,5)A with a resolution of 0.2 A
        # this should be left to the user eventually
        grid_xyz_coords = []
        for i, j, k in product(
            range(resolution),
            range(resolution),
            range(resolution),
        ):
            itrans = -5 + 0.2 * i
            jtrans = -5 + 0.2 * j
            ktrans = -5 + 0.2 * k
            grid_xyz_coords.append([itrans, jtrans, ktrans])
        # write the grid to a .dat file
        with open(output_directory.joinpath("grid.dat"), "w") as file:
            for xyz in grid_xyz_coords:
                for c in xyz:
                    file.write(str(c) + " ")
                file.write("\n")

    def calculate_electrostatic_potential(
        self,
        output_directory: pathlib.Path | str,
        resolution: int = 51,
        num_threads: int = 14,
        optimize: bool = False,
    ) -> None:

        output_directory = pathlib.Path(output_directory)
        _create_directory(output_directory)

        self._generate_voxel_grid(resolution, output_directory)
        original_directory = os.getcwd()  # aka OGD
        os.chdir(output_directory)

        psi4.set_options(
            {
                "basis": "aug-cc-pVDZ",
                "CUBEPROP_TASKS": ["ESP"],
                "CUBEPROP_FILEPATH": str(output_directory),
                "reference": "uhf",
            }
        )
        psi4.core.set_num_threads(num_threads)

        psi4_mol = psi4.core.Molecule.from_arrays(
            self._coordinates, elem=self._elements
        )
        psi4.core.set_output_file(
            str(output_directory.joinpath("output.dat")), False
        )
        self.output = output_directory.joinpath("output.dat")

        if optimize:
            print("optimizing!")
            psi4.optimize("PBE", molecule=psi4_mol)

        print("calculating ESP")
        E, wfn = psi4.prop(
            "PBE", molecule=psi4_mol, properties=["GRID_ESP"], return_wfn=True
        )
        psi4.cubeprop(wfn)

        os.chdir(original_directory)
