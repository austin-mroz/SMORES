import pathlib

import smores
import smores.psi4


def test_get_steric_parameters(tmp_path: pathlib.Path) -> None:
    molecule = smores.rdkit_from_smiles("Br")
    smores.psi4.calculate_electrostatic_potential(
        molecule=molecule,
        output_directory=tmp_path,
        grid_origin=(-5, -5, -5),
        grid_length=10,
        num_voxels_per_dimension=50,
    )
    esp_molecule = smores.EspMolecule.from_cube_file(
        path=tmp_path / "ESP.cube",
    )
    steric_parameters = esp_molecule.get_steric_parameters(0, 1)
