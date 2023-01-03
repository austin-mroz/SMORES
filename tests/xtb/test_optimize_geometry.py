import pathlib

import numpy as np
import rdkit.Chem.AllChem as rdkit

import smores


def test_optimize_geometry(tmp_path: pathlib.Path) -> None:
    molecule = smores.rdkit_from_smiles("CBr")
    optimized = smores.xtb.optimize_geometry(molecule, tmp_path / "xtb")
    assert rdkit.MolToSmiles(molecule) == rdkit.MolToSmiles(optimized)
    assert not np.all(
        np.isclose(
            molecule.GetConformer().GetPositions(),
            optimized.GetConformer().GetPositions(),
        ),
    )
