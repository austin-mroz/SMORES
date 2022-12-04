import pathlib
import subprocess

import rdkit.Chem.AllChem as rdkit

from smores._internal.utilities import current_working_directory


def optimize_geometry(
    molecule: rdkit.Mol,
    output_directory: pathlib.Path | str,
    level: str = "extreme",
) -> pathlib.Path:
    """
    Optimize the geometry of a molecule.

    """
    output_directory = pathlib.Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    xyz_path = str(output_directory / "geom.xyz")
    rdkit.MolToXYZFile(molecule, xyz_path)

    with current_working_directory(output_directory):
        subprocess.run(
            [
                "xtb",
                xyz_path,
                "--opt",
                level,
            ],
            check=True,
        )

    return rdkit.MolFromXYZFile(str(output_directory / "xtbopt.xyz"))
