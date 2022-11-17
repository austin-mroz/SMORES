import pathlib
import typing
import rdkit.Chem.AllChem as rdkit
import subprocess


def optimize_geometry(
        molecule: rdkit.Mol,
        output_directory: pathlib.Path | str,
        level: str = "extreme",
) -> pathlib.Path:
    """
    Calculate the electrostatic potential...

    """
    output_directory = pathlib.Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    xyz_path = str(output_directory / "geom.xyz")
    rdkit.MolToXYZFile(molecule, xyz_path)

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
