import pathlib
import subprocess

import rdkit.Chem.AllChem as rdkit

from smores._internal.utilities import current_working_directory


def optimize_geometry(
    molecule: rdkit.Mol,
    output_directory: pathlib.Path | str,
    level: str = "extreme",
    xtb_path: pathlib.Path | str = "xtb",
) -> rdkit.Mol:
    """
    Optimize the geometry of a molecule.

    .. warning::

        xtb will not work if you have an activate Python environment
        in which ``psi4`` is installed.

    Examples:

        *Optimize the geometry of a molecule*


        .. testcode:: optimize-the-geometry-of-a-molecule
            :hide:

            import os
            import tempfile
            tmp_dir = tempfile.TemporaryDirectory()
            os.chdir(tmp_dir.name)

        .. testcode:: optimize-the-geometry-of-a-molecule

            import smores
            molecule = smores.rdkit_from_smiles("CBr")
            optimized = smores.xtb.optimize_geometry(molecule, "CBr_xtb")

    Parameters:
        molecule:
            The molecule to optimize.
        output_directory:
            The directory in which xtb places its files.
        level:
            The optimization level. Passed directly to xtb.
        xtb_path:
            The path to the xtb binary.
    Returns:
        The optimized molecule.

    """
    output_directory = pathlib.Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    xyz_path = str(output_directory / "geom.xyz")
    rdkit.MolToXYZFile(molecule, xyz_path)

    with current_working_directory(output_directory):
        process = subprocess.run(
            [
                str(xtb_path),
                xyz_path,
                "--opt",
                level,
            ],
            check=True,
            capture_output=True,
            text=True,
        )

    with open(output_directory / "xtb.stdout", "w") as f:
        f.write(process.stdout)

    with open(output_directory / "xtb.stderr", "w") as f:
        f.write(process.stderr)

    xtb_molecule = rdkit.MolFromXYZFile(str(output_directory / "xtbopt.xyz"))
    optimized = rdkit.Mol(molecule)
    optimized.RemoveConformer(0)
    optimized.AddConformer(xtb_molecule.GetConformer())
    return optimized
