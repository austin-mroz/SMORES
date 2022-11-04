from __future__ import annotations

import numpy as np
import numpy.typing as mpt
from openbabel import openbabel
import pathlib
import rdkit.Chem.AllChem as rdkit
from itertools import product
import pandas as pd
import morfeus
import dbstep.Dbstep as db
from dataclasses import dataclass

from . import utilities


@dataclass(frozen=True)
class SterimolParameters:
    L: float
    B1: float
    B5: float


class SterimolCalculator:
    """
    handles calculation of sterimol parameters using
    variety of methods and softwares
    """
    
    def __init__(
            self,
    ) -> None:
        self._molecule_xyz_path = None
        self._molecule_mol_path = None


    @classmethod
    def init_from_file(
            cls,
            molecule_path: pathlib.Path | str,
    ) -> SterimolCalculator:

        instance = cls.__new__(cls)

        if isinstance(molecule_path, str):
            molecule_path = pathlib.Path(molecule_path)

        if str(molecule_path.parts[-1])[-3:] == 'xyz':
            instance._molecule_mol_path = _convert_molecule_file_type(
                    molecule_path,
                    conv_from='xyz',
                    conv_to='mol',
            )
            instance._molecule_xyz_path = molecule_path
        elif str(molecule_path.parts[-1])[-3:] == 'mol':
            instance._molecule_mol_path = molecule_path
            instance._molecule_xyz_path = _convert_molecule_file_type(
                    molecule_path,
                    conv_from='mol',
                    conv_to='xyz',
            )
        else:
            raise ValueError(f'{molecule_path.parts[-1]} is not an accepted format. Please use mol or xyz file type')

        return instance


    def set_core_smiles(
            self,
            core_smiles: str,
    ) -> None:
        self.core_smiles = core_smiles


    def set_atom_idx(
            self,
            atom_1_idx: int,
            atom_2_idx: int,
    ) -> None:

        self.connecting_atom_idx_dict = {
                'core': atom_1_idx,
                'substituent': atom_2_idx,
        }


    def calculate_sterimol_parameters(
            self,
            method_types: list[str],
            radii_types: list[str] = ['streusel'],
    ) -> None:

        if hasattr(self, "core_smiles"):
            self.connecting_atom_idx = _get_sterimol_atom_idx(self)
        if not hasattr(self, "connecting_atom_idx_dict"):
            raise ValueError("You must set either the atom idx or core smiles.")

        accepted_method_types = set(['db','morfeus'])
        sterimol_parameters = []

        method_types = [method_type.lower() for method_type in method_types]

        if any(method_type in method_types for method_type in accepted_method_types):
            if 'morfeus' in method_types:
                sterimol_parameters.append(
                        _collate_morfeus_sterimol_parameters(
                            self._molecule_xyz_path,
                            radii_types,
                            self.connecting_atom_idx_dict,
                        )
                )
            if 'db' in method_types:
                sterimol_parameters.append(
                        _calculate_dbstep_parameters(
                            self._molecule_xyz_path,
                            self.connecting_atom_idx_dict,
                        )
                )
        else:
            raise NameError(f'We only support {accepted_method_types}')

        return sterimol_parameters


@dataclass(frozen=True)
class SterimolParameter:
    method: str
    radii_type: str
    L: float
    B1: float
    B5: float


def _calculate_dbstep_parameters(
        molecule_xyz_path: pathlib.Path,
        connecting_atom_idx: dict,
) -> list[SterimolParameters]:
    dbstep_sterimol = db.dbstep(
            str(molecule_xyz_path),
            atom1=int(connecting_atom_idx['core']),
            atom2=int(connecting_atom_idx['substituent']),
            sterimol=True,
            measure='classic',
    )
    return [SterimolParameter(
            method='DBStep',
            radii_type='vdW',
            L=dbstep_sterimol.L,
            B1=dbstep_sterimol.Bmin,
            B5=dbstep_sterimol.Bmax,
    )]


def _collate_morfeus_sterimol_parameters(
        molecule_xyz_path: pathlib.Path,
        radii_types: list[str],
        connecting_atom_idx: dict,
) -> list[SterimolParameters]:
    parameters = []
    for radii_type in radii_types:
        morfeus_sterimol = _calculate_morfeus_sterimol_parameters(
                molecule_xyz_path,
                radii_type,
                connecting_atom_idx['core'],
                connecting_atom_idx['substituent'],
        )
        parameters.append(
                SterimolParameter(
                    method='morfeus',
                    radii_type=radii_type,
                    L=morfeus_sterimol.L_value,
                    B1=morfeus_sterimol.B_1_value,
                    B5=morfeus_sterimol.B_5_value,
                )
        )
    return parameters


def _calculate_morfeus_sterimol_parameters(
        molecule_xyz_path: pathlib.Path,
        radii_type: str,
        core_atom_idx: int,
        subst_atom_idx: int,
) -> morfeus.Sterimol:

    xyz_data = utilities.read_xyz(molecule_xyz_path)
    if radii_type == 'streusel':
        streusel_radii = utilities.get_streusel_radii(xyz_data)
        morfeus_sterimol = morfeus.Sterimol(
                xyz_data.elements,
                xyz_data.coordinates,
                core_atom_idx,
                subst_atom_idx,
                radii=streusel_radii,
        )
    else:
        morfeus_sterimol = morfeus.Sterimol(
                xyz_data.elements,
                xyz_data.coordinates,
                core_atom_idx,
                subst_atom_idx,
                radii_type=radii_type,
        )
    return morfeus_sterimol


def _get_sterimol_atom_idx(
        self,
) -> dict:
    mol = rdkit.MolFromMolFile(str(self._molecule_mol_path))

    if not hasattr(self, 'core_smiles'):
        raise NameError('core smiles has not been defined. Set core_smiles using > SterimolCalculator().set_core_smiles(core_smiles).')
    core_pattern = rdkit.MolFromSmiles(self.core_smiles)
    core_matches = mol.GetSubstructMatch(core_pattern)

    core_coords = []
    core_elements = []
    for core_match in core_matches:
        core_elements.append(mol.GetAtoms()[core_match].GetSymbol())
        x, y, z = list(mol.GetConformer().GetAtomPosition(core_match))
        core_coords.append([x, y, z])
    core_xyz = utilities.XyzData(
            elements=core_elements,
            coordinates=core_coords,
    )
    mol_xyz = utilities.read_xyz(self._molecule_xyz_path)
    print(mol_xyz)
    print(core_xyz)
    exit()
    core_xyz_array = utilities.get_full_xyz_array(core_xyz)
    mol_xyz_array = utilities.get_full_xyz_array(mol_xyz, removeHs=True)

    mol_C_coords_array = np.round(np.delete(mol_xyz_array, 0, axis=1).astype(np.float32),3)
    core_C_coords_array = np.round(np.delete(core_xyz_array, 0, axis=1).astype(np.float32),3)

    mol_C_coords_df = pd.DataFrame(mol_C_coords_array)
    core_C_coords_df = pd.DataFrame(core_C_coords_array)

    # the difference between the full molecule (mol) and the core
    # is the substituent
    subst_coords_df = pd.concat([mol_C_coords_df, core_C_coords_df]).drop_duplicates(keep=False)
    subst_coords = np.asarray(subst_coords_df)

    # we want the core Carbon and substituent Carbon that are closest together
    for core_coord, subst_coord in product(core_C_coords_array, subst_coords):
        if np.linalg.norm(core_coord.reshape(-1,1).T-subst_coord) < 1.6:
            core_min_C_coord = core_coord
            subst_min_C_coord = subst_coord
            break
    core_connecting_atom_idx = np.where(np.all(core_min_C_coord == np.round(mol_xyz.coordinates,3),axis=1))
    subst_connecting_atom_idx = np.where(np.all(subst_min_C_coord == np.round(mol_xyz.coordinates,3),axis=1))

    connecting_atom_idx_dict = {
            'core': core_connecting_atom_idx[0],
            'substituent': subst_connecting_atom_idx[0]
    }

    return connecting_atom_idx_dict


def _convert_molecule_file_type(
        input_path: pathlib.Path,
        conv_from: str,
        conv_to: str,
) -> pathlib.Path:
    calculation_directory = input_path.parents[0]
    input_name = input_path.parts[-1]
    output_name = input_name.removesuffix(conv_from) + conv_to

    output_path = calculation_directory.joinpath(output_name)

    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats(conv_from, conv_to)

    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, str(input_path))
    obConversion.WriteFile(mol, str(output_path))

    return output_path
