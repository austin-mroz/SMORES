from __future__ import annotations

import os
import pathlib
from itertools import product
import logging

import numpy as np
import numpy.typing as npt
import psi4
import rdkit.Chem.AllChem as AllChem

from . import utilities
from . import constants


class InvalidDirectoryError(Exception):
    pass


class Molecule:
    """
    class to handle molecule objects for featurization
    """

    def __init__(
        self,
        smiles: str,
    ) -> None:
        rdkit_mol = AllChem.AddHs(AllChem.MolFromSmiles(smiles))
        self._elements = [atom.GetSymbol() for atom in rdkit_mol.GetAtoms()]
        AllChem.EmbedMolecule(rdkit_mol)
        rdkit_coordinates = (
            rdkit_mol.GetConformer(0).GetPositions().astype(np.float32)
        )
        center_of_mass = _get_center_of_mass(self._elements, rdkit_coordinates)
        self._coordinates = rdkit_coordinates - center_of_mass


    @classmethod
    def init_from_xyz(
        cls,
        path: pathlib.Path | str,
    ) -> Molecule:

        if isinstance(path, str):
            path = pathlib.Path(path)

        # Create an instance without calling the __init__ method.
        instance = cls.__new__(cls)
        xyz_data = utilities.read_xyz(path)
        instance._elements = xyz_data.elements
        center_of_mass = _get_center_of_mass(
            instance._elements,
            xyz_data.coordinates,
        )
        instance._coordinates = xyz_data.coordinates - center_of_mass
        return instance


    def generate_xtb_starting_structure(
            self,
            output_directory: pathlib.Path | str,
            xyz_file_name: str = "initial_structure",
    ) -> None:

        # create the output directory
        if isinstance(output_directory, str):
            output_directory = pathlib.Path(output_directory)
        _create_directory(output_directory)

        # write the xyz file with the updated coordinates
        xyz_data = utilities.XyzData(
                elements=self._elements,
                coordinates=self._coordinates,
        )
        xyz_path = output_directory.joinpath(f"{xyz_file_name}.xyz")
        utilities.write_xyz(
                xyz_path,
                xyz_data,
        )


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


def _create_directory(path: pathlib.Path) -> None:
    if path.exists() and not path.is_dir():
        raise InvalidDirectoryError(f"{path} is not a valid directory.")

    if not path.exists():
        path.mkdir(exist_ok=True, parents=True)


def _get_center_of_mass(
    elements: list[str],
    coordinates: npt.NDArray[np.float32],
) -> npt.NDArray[np.float32]:
    atom_masses = np.array(
        [constants.atomic_mass[element] for element in elements]
    )
    scaled_coordinates = coordinates * atom_masses[:, np.newaxis]
    return scaled_coordinates.sum(axis=0) / atom_masses.sum()
