import pathlib
from dataclasses import dataclass

import morfeus
import numpy as np
import numpy.typing as npt
import pytest
import rdkit.Chem.AllChem as rdkit
from pytest_lazyfixture import lazy_fixture

import smores


@dataclass(frozen=True, slots=True)
class CaseData:
    atoms: tuple[str, ...]
    positions: npt.NDArray[np.float32]
    molecule: smores.Molecule


def test_smores_parameters_match_sterimol_if_same_radii_are_used(
    case_data: CaseData,
) -> None:

    dummy_index = 0
    attached_index = 1
    radii = {"H": 1.0, "Br": 1.0}
    sterimol = morfeus.Sterimol(
        elements=case_data.atoms,
        coordinates=case_data.positions,
        dummy_index=dummy_index + 1,
        attached_index=attached_index + 1,
        radii=tuple(radii.values()),
    )
    params = case_data.molecule.get_steric_parameters(
        dummy_index=dummy_index,
        attached_index=attached_index,
        radii=radii,
    )
    assert params.L == sterimol.L_value
    assert params.B1 == sterimol.B_1_value
    assert params.B5 == sterimol.B_5_value


@pytest.fixture(
    params=(
        lazy_fixture("default_init_molecule"),
        lazy_fixture("molecule_from_xyz_file"),
        lazy_fixture("molecule_from_mol_file"),
        lazy_fixture("molecule_from_smiles"),
    ),
)
def case_data(request) -> CaseData:
    return request.param


@pytest.fixture
def default_init_molecule(rdkit_molecule: rdkit.Mol) -> CaseData:
    return CaseData(
        atoms=tuple(atom.GetSymbol() for atom in rdkit_molecule.GetAtoms()),
        positions=rdkit_molecule.GetConformer(0).GetPositions(),
        molecule=smores.Molecule(
            atoms=tuple(
                atom.GetSymbol() for atom in rdkit_molecule.GetAtoms()
            ),
            positions=rdkit_molecule.GetConformer(0).GetPositions(),
        ),
    )


@pytest.fixture
def molecule_from_xyz_file(
    tmp_path: pathlib.Path,
    rdkit_molecule: rdkit.Mol,
) -> CaseData:

    rdkit.MolToXYZFile(rdkit_molecule, str(tmp_path / "molecule.xyz"))
    return CaseData(
        atoms=tuple(atom.GetSymbol() for atom in rdkit_molecule.GetAtoms()),
        positions=rdkit_molecule.GetConformer(0).GetPositions(),
        molecule=smores.Molecule.from_xyz_file(tmp_path),
    )


@pytest.fixture
def molecule_from_mol_file(
    tmp_path: pathlib.Path,
    rdkit_molecule: rdkit.Mol,
) -> CaseData:

    rdkit.MolToMolFile(rdkit_molecule, str(tmp_path / "molecule.mol"))
    return CaseData(
        atoms=tuple(atom.GetSymbol() for atom in rdkit_molecule.GetAtoms()),
        positions=rdkit_molecule.GetConformer(0).GetPositions(),
        molecule=smores.Molecule.from_mol_file(tmp_path),
    )


@pytest.fixture
def molecule_from_smiles(rdkit_molecule: rdkit.Mol) -> CaseData:
    return CaseData(
        atoms=tuple(atom.GetSymbol() for atom in rdkit_molecule.GetAtoms()),
        positions=rdkit_molecule.GetConformer(0).GetPositions(),
        molecule=smores.Molecule.from_smiles(
            smiles=rdkit.MolToSmiles(rdkit_molecule),
            positions=rdkit_molecule.GetConformer(0).GetPositions(),
        ),
    )


@pytest.fixture
def rdkit_molecule() -> rdkit.Mol:
    molecule = rdkit.MolFromSmiles("Br")
    rdkit.EmbedMolecule(molecule, rdkit.ETKDGv3())
    return molecule
