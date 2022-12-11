import pathlib

import smores
import smores.psi4


def test_get_steric_parameters(datadir: pathlib.Path) -> None:
    esp_molecule = smores.EspMolecule.from_cube_file(path=datadir / "ESP.cube")
    steric_parameters = esp_molecule.get_steric_parameters(0, 1)
    print(steric_parameters)
