import pathlib
import subprocess

import rdkit.Chem.AllChem as rdkit

from smores._internal.utilities import current_working_directory


def optimize_geometry(
    molecule: rdkit.Mol,
    output_directory: pathlib.Path | str,
    level: str = "extreme",
) -> rdkit.Mol:
    """
    Optimize the geometry of a molecule.

    Examples:
        *Optimize the geometry of a molecule*

    Parameters:
        molecule:
            The molecule to optimize.
        output_directory:
            The directory in which xtb places its files.
        level:
            The optimization level. Passed directly to xtb.
    Returns:
        The optimized molecule.

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
