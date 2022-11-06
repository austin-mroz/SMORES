import pytest
import rdkit.Chem.AllChem as rdkit

import smores


@pytest.mark.parametrize(
    ("core_smiles", "substituent_smiles", "expected_product_smiles"),
    (
        ("CC[Sn]Br", "CC[Ge]Br", "CC[Sn][Ge]CC"),
        ("Br[Sn]CC", "Br[Ge]CC", "CC[Sn][Ge]CC"),
        ("CC[Sn](Br)CC", "CC[Ge](Br)CC", "CC[Sn](CC)[Ge](CC)CC"),
    ),
)
def test_combine(
    core_smiles: str,
    substituent_smiles: str,
    expected_product_smiles: str,
) -> None:

    core = smores.rdkit_from_smiles(core_smiles)
    substituent = smores.rdkit_from_smiles(substituent_smiles)
    (combo,) = smores.combine([core], [substituent])
    expected_dummy_index = next(
        atom.GetIdx()
        for atom in combo.product.GetAtoms()
        if atom.GetSymbol() == "Sn"
    )
    expected_attached_index = next(
        atom.GetIdx()
        for atom in combo.product.GetAtoms()
        if atom.GetSymbol() == "Ge"
    )
    assert combo.core is core
    assert combo.substituent is substituent

    rdkit.SanitizeMol(combo.product)
    expected_product_smiles = _to_canonical_smiles(expected_product_smiles)
    product_smiles = rdkit.MolToSmiles(
        rdkit.RemoveHs(combo.product),
        canonical=True,
    )

    assert product_smiles == expected_product_smiles
    assert combo.dummy_index == expected_dummy_index
    assert combo.attached_index == expected_attached_index


def _to_canonical_smiles(smiles: str) -> str:
    return rdkit.MolToSmiles(rdkit.MolFromSmiles(smiles), canonical=True)
