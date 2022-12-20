import pathlib
import typing
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
    radii: tuple[float, ...]
    molecule: smores.Molecule


def test_smores_parameters_match_sterimol_if_same_radii_are_used(
    case_data: CaseData,
) -> None:

    sterimol = morfeus.Sterimol(
        elements=case_data.atoms,
        coordinates=case_data.positions,
        dummy_index=case_data.molecule.get_dummy_index() + 1,
        attached_index=case_data.molecule.get_attached_index() + 1,
        radii=case_data.radii,
    )
    params = case_data.molecule.get_steric_parameters()
    assert params.L == pytest.approx(sterimol.L_value, abs=1e-4)
    assert params.B1 == pytest.approx(sterimol.B_1_value)
    assert params.B5 == pytest.approx(sterimol.B_5_value)


@pytest.fixture(
    params=(
        lazy_fixture("default_init_molecule"),
        lazy_fixture("molecule_from_xyz_file"),
        lazy_fixture("molecule_from_mol_file"),
        lazy_fixture("molecule_from_smiles"),
    ),
)
def case_data(request: typing.Any) -> CaseData:
    return request.param


@pytest.fixture
def default_init_molecule(rdkit_molecule: rdkit.Mol) -> CaseData:
    radii = (1.0, 1.0)
    return CaseData(
        atoms=tuple(atom.GetSymbol() for atom in rdkit_molecule.GetAtoms()),
        positions=rdkit_molecule.GetConformer(0).GetPositions(),
        radii=radii,
        molecule=smores.Molecule(
            atoms=tuple(
                atom.GetSymbol() for atom in rdkit_molecule.GetAtoms()
            ),
            dummy_index=0,
            attached_index=1,
            positions=rdkit_molecule.GetConformer(0).GetPositions(),
            radii=radii,
        ),
    )


@pytest.fixture
def molecule_from_xyz_file(
    tmp_path: pathlib.Path,
    rdkit_molecule: rdkit.Mol,
) -> CaseData:

    xyz_file = tmp_path / "molecule.xyz"
    radii = (1.0, 1.0)
    rdkit.MolToXYZFile(rdkit_molecule, str(xyz_file))
    return CaseData(
        atoms=tuple(atom.GetSymbol() for atom in rdkit_molecule.GetAtoms()),
        positions=rdkit_molecule.GetConformer(0).GetPositions(),
        radii=radii,
        molecule=smores.Molecule.from_xyz_file(
            path=xyz_file,
            dummy_index=0,
            attached_index=1,
            radii=radii,
        ),
    )


@pytest.fixture
def molecule_from_mol_file(
    tmp_path: pathlib.Path,
    rdkit_molecule: rdkit.Mol,
) -> CaseData:

    mol_file = tmp_path / "molecule.mol"
    radii = (1.0, 1.0)
    rdkit.MolToMolFile(rdkit_molecule, str(mol_file))
    return CaseData(
        atoms=tuple(atom.GetSymbol() for atom in rdkit_molecule.GetAtoms()),
        positions=rdkit_molecule.GetConformer(0).GetPositions(),
        radii=radii,
        molecule=smores.Molecule.from_mol_file(
            path=mol_file,
            dummy_index=0,
            attached_index=1,
            radii=radii,
        ),
    )


@pytest.fixture
def molecule_from_smiles(rdkit_molecule: rdkit.Mol) -> CaseData:
    radii = (1.0, 1.0)
    return CaseData(
        atoms=tuple(atom.GetSymbol() for atom in rdkit_molecule.GetAtoms()),
        positions=rdkit_molecule.GetConformer(0).GetPositions(),
        radii=radii,
        molecule=smores.Molecule.from_smiles(
            smiles=rdkit.MolToSmiles(rdkit_molecule),
            dummy_index=0,
            attached_index=1,
            positions=rdkit_molecule.GetConformer(0).GetPositions(),
            radii=radii,
        ),
    )


@pytest.fixture
def rdkit_molecule() -> rdkit.Mol:
    molecule = rdkit.AddHs(rdkit.MolFromSmiles("Br"))
    rdkit.EmbedMolecule(molecule, rdkit.ETKDGv3())
    return molecule
