import pytest
import rdkit.Chem.AllChem as rdkit

import smores


def test_combine(
    core_smiles: str,
    substituent_smiles: str,
) -> None:

    core = smores.rdkit_from_smiles(core_smiles)
    substituent = smores.rdkit_from_smiles(substituent_smiles)
    (combo,) = smores.combine([core], [substituent])
    expected_dummy_index = next(
        atom.GetIdx()
        for atom in combo.product.GetAtoms()
        if atom.GetSymbol() == "Si"
    )
    expected_attached_index = next(
        atom.GetIdx()
        for atom in combo.product.GetAtoms()
        if atom.GetSymbol() == "Ge"
    )
    assert core is combo.core
    assert substituent is combo.substituent

    expected_product_smiles = _get_expected_product_smiles(
        core_smiles, substituent_smiles
    )
    assert expected_product_smiles == rdkit.MolToSmiles(
        combo.product,
        canonical=True,
    )

    assert expected_dummy_index == combo.dummy_index
    assert expected_attached_index == combo.attached_index


def _get_expected_product_smiles(
    core_smiles: str,
    substituent_smiles: str,
) -> str:

    product_smiles = core_smiles.replace(
        "Br",
        substituent_smiles.replace("Br", "").replace("()", ""),
    )
    return rdkit.MolToSmiles(
        rdkit.MolFromSmiles(product_smiles),
        canonical=True,
    )


@pytest.fixture(
    params=(
        "CCC[Si]Br",
        "Br[Si]CCC",
        "CCC[Si](Br)CCC",
        # "c1cccc[si]1Br"
    ),
)
def core_smiles(request) -> str:
    return request.param


@pytest.fixture(
    params=(
        "CCC[Ge]Br",
        "Br[Ge]CCC",
        "CCC[Ge](Br)CCC",
    ),
)
def substituent_smiles(request) -> str:
    return request.param
