def _generate_voxel_grid(
    self,
    resolution: int,
    output_directory: pathlib.Path,
) -> None:
    # we want to make a box with one corner at (-5,5,5)A and one at
    # (5,5,5)A with a resolution of 0.2 A
    # this should be left to the user eventually
    grid_xyz_coords = []
    for i, j, k in product(
        range(resolution),
        range(resolution),
        range(resolution),
    ):
        itrans = -5 + 0.2 * i
        jtrans = -5 + 0.2 * j
        ktrans = -5 + 0.2 * k
        grid_xyz_coords.append([itrans, jtrans, ktrans])
    # write the grid to a .dat file
    with open(output_directory.joinpath("grid.dat"), "w") as file:
        for xyz in grid_xyz_coords:
            for c in xyz:
                file.write(str(c) + " ")
            file.write("\n")


def calculate_electrostatic_potential(
    self,
    output_directory: pathlib.Path | str,
    resolution: int = 51,
    num_threads: int = 14,
    optimize: bool = False,
) -> None:

    output_directory = pathlib.Path(output_directory)
    _create_directory(output_directory)

    self._generate_voxel_grid(resolution, output_directory)
    original_directory = os.getcwd()  # aka OGD
    os.chdir(output_directory)

    psi4.set_options(
        {
            "basis": "aug-cc-pVDZ",
            "CUBEPROP_TASKS": ["ESP"],
            "CUBEPROP_FILEPATH": str(output_directory),
            "reference": "uhf",
        }
    )
    psi4.core.set_num_threads(num_threads)

    psi4_mol = psi4.core.Molecule.from_arrays(
        self._coordinates, elem=self._elements
    )
    psi4.core.set_output_file(
        str(output_directory.joinpath("output.dat")), False
    )
    self.output = output_directory.joinpath("output.dat")

    if optimize:
        print("optimizing!")
        psi4.optimize("PBE", molecule=psi4_mol)

    print("calculating ESP")
    E, wfn = psi4.prop(
        "PBE", molecule=psi4_mol, properties=["GRID_ESP"], return_wfn=True
    )
    psi4.cubeprop(wfn)

    os.chdir(original_directory)
