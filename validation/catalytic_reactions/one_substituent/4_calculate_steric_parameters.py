#!python
import argparse
import csv
import json
import pathlib
import typing
from dataclasses import dataclass

import flour
import morfeus
import numpy as np
import numpy.typing as npt
import pyvista as pv
import scipy

import smores

_OUTPUT_CSV_COLUMNS = (
    "name",
    "core",
    "substituent",
    "smiles",
    "xyz_file",
    "fragments_file",
    "dummy_index",
    "attached_index",
    "radii_type",
    "L",
    "B1",
    "B5",
)


def main() -> None:
    args = _get_command_line_arguments()
    args.output_directory.mkdir(parents=True, exist_ok=True)

    radii_types = ("alvarez", "bondi", "crc", "rahm", "pyykko", "truhlar")
    with open(
        args.output_directory / "steric_parameters.csv", "w"
    ) as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=_OUTPUT_CSV_COLUMNS)
        writer.writeheader()

        for catalyst_input_file in args.input_file:
            for row in tuple(_get_rows(catalyst_input_file)):
                print(row.xyz_file)
                smores_molecule = smores.Molecule.from_xyz_file(
                    path=row.xyz_file,
                    dummy_index=row.dummy_index,
                    attached_index=row.attached_index,
                )
                smores_params = smores_molecule.get_steric_parameters()
                writer.writerow(
                    {
                        "name": row.name,
                        "core": row.core,
                        "substituent": row.substituent,
                        "smiles": row.smiles,
                        "xyz_file": row.xyz_file,
                        "fragments_file": row.fragments_file,
                        "dummy_index": row.dummy_index,
                        "attached_index": row.attached_index,
                        "radii_type": "streusel",
                        "L": smores_params.L,
                        "B1": smores_params.B1,
                        "B5": smores_params.B5,
                    }
                )

                with open(row.fragments_file) as f:
                    indices_dict = json.load(f)
                    core_indices = indices_dict["core_indices"]
                    substituent_indices = indices_dict["substituent_indices"]

                smores_core_excluded_molecule = smores.Molecule.from_xyz_file(
                    path=row.xyz_file,
                    dummy_index=row.dummy_index,
                    attached_index=row.attached_index,
                    excluded_indices=core_indices,
                )
                smores_core_excluded_params = (
                    smores_core_excluded_molecule.get_steric_parameters()
                )
                writer.writerow(
                    {
                        "name": row.name,
                        "core": row.core,
                        "substituent": row.substituent,
                        "smiles": row.smiles,
                        "xyz_file": row.xyz_file,
                        "fragments_file": row.fragments_file,
                        "dummy_index": row.dummy_index,
                        "attached_index": row.attached_index,
                        "radii_type": "streusel_core_excluded",
                        "L": smores_core_excluded_params.L,
                        "B1": smores_core_excluded_params.B1,
                        "B5": smores_core_excluded_params.B5,
                    }
                )
                """
                smores_esp_molecule = smores.EspMolecule.from_cube_file(
                    path=row.xyz_file.parent / "ESP.cube",
                    dummy_index=row.dummy_index,
                    attached_index=row.attached_index,
                    included_indices=substituent_indices,
                )

                electric_field_surface = (
                    smores_esp_molecule.get_electric_field_surface()
                )
                flour.write_cube(
                    path=row.xyz_file.parent / "STREUSEL_core_excluded.cube",
                    title1="t1",
                    title2="t2",
                    atoms=smores_esp_molecule._atoms,
                    charges=np.zeros(
                        len(smores_esp_molecule._atoms), dtype=np.float64
                    ),
                    positions=smores_esp_molecule._positions,
                    voxel_origin=electric_field_surface.voxel_origin,
                    voxel_size=np.identity(3)
                    * electric_field_surface.voxel_size,
                    voxels=np.array(
                        electric_field_surface.voxels, dtype=np.float64
                    ),
                )

                _plot_surface(
                    row.xyz_file.parent / "STREUSEL_core_excluded.cube",
                    attached_atom_idx=row.attached_index,
                    dummy_atom_idx=row.dummy_index,
                )

                esp_smores_params = smores_esp_molecule.get_steric_parameters()
                writer.writerow(
                    {
                        "name": row.name,
                        "core": row.core,
                        "substituent": row.substituent,
                        "smiles": row.smiles,
                        "xyz_file": row.xyz_file,
                        "fragments_file": row.fragments_file,
                        "dummy_index": row.dummy_index,
                        "attached_index": row.attached_index,
                        "radii_type": "streusel_cube_core_excluded",
                        "L": esp_smores_params.L,
                        "B1": esp_smores_params.B1,
                        "B5": esp_smores_params.B5,
                    }
                )
                """
                smores_esp_molecule = smores.EspMolecule.from_cube_file(
                    path=row.xyz_file.parent / "ESP.cube",
                    dummy_index=row.dummy_index,
                    attached_index=row.attached_index,
                )
                esp_smores_params = smores_esp_molecule.get_steric_parameters(
                    plot=True,
                    output_path=args.output_directory
                    / f"{row.name}_pointcloud.png",
                )
                writer.writerow(
                    {
                        "name": row.name,
                        "core": row.core,
                        "substituent": row.substituent,
                        "smiles": row.smiles,
                        "xyz_file": row.xyz_file,
                        "fragments_file": row.fragments_file,
                        "dummy_index": row.dummy_index,
                        "attached_index": row.attached_index,
                        "radii_type": "streusel_cube_core_included",
                        "L": esp_smores_params.L * 0.529,
                        "B1": esp_smores_params.B1 * 0.529,
                        "B5": esp_smores_params.B5 * 0.529,
                    }
                )

                for radii_type in radii_types:
                    sterimol = _get_sterimol(row, radii_type)
                    writer.writerow(
                        {
                            "name": row.name,
                            "core": row.core,
                            "substituent": row.substituent,
                            "smiles": row.smiles,
                            "xyz_file": row.xyz_file,
                            "fragments_file": row.fragments_file,
                            "dummy_index": row.dummy_index,
                            "attached_index": row.attached_index,
                            "radii_type": radii_type,
                            "L": sterimol.L_value,
                            "B1": sterimol.B_1_value,
                            "B5": sterimol.B_5_value,
                        },
                    )
                    sterimol_excluded = _get_sterimol(
                        molecule=row,
                        radii_type=radii_type,
                        excluded_indices=core_indices,
                    )
                    writer.writerow(
                        {
                            "name": row.name,
                            "core": row.core,
                            "substituent": row.substituent,
                            "smiles": row.smiles,
                            "xyz_file": row.xyz_file,
                            "fragments_file": row.fragments_file,
                            "dummy_index": row.dummy_index,
                            "attached_index": row.attached_index,
                            "radii_type": f"{radii_type}_core_excluded",
                            "L": sterimol_excluded.L_value,
                            "B1": sterimol_excluded.B_1_value,
                            "B5": sterimol_excluded.B_5_value,
                        },
                    )


@dataclass(frozen=True, slots=True)
class CsvRow:
    name: str
    core: str
    substituent: str
    smiles: str
    xyz_file: pathlib.Path
    fragments_file: pathlib.Path
    dummy_index: int
    attached_index: int


def _get_rows(
    path: pathlib.Path,
) -> typing.Iterator[CsvRow]:
    with open(path) as csv_file:
        reader = csv.DictReader(csv_file)
        for row in reader:
            yield CsvRow(
                name=row["name"],
                core=row["core"],
                substituent=row["substituent"],
                smiles=row["smiles"],
                xyz_file=pathlib.Path(row["xyz_file"]),
                fragments_file=pathlib.Path(row["fragments_file"]),
                dummy_index=int(row["dummy_index"]),
                attached_index=int(row["attached_index"]),
            )


def _get_sterimol(
    molecule: CsvRow,
    radii_type: str,
    excluded_indices: list[int] | None = None,
) -> morfeus.Sterimol:
    elements, coordinates = morfeus.read_xyz(molecule.xyz_file)
    return morfeus.Sterimol(
        elements=elements,
        coordinates=coordinates,
        dummy_index=molecule.dummy_index + 1,
        attached_index=molecule.attached_index + 1,
        radii_type=radii_type,
        excluded_atoms=(
            None
            if excluded_indices is None
            else [index + 1 for index in excluded_indices]
        ),
    )


def _plot_surface(
    cube_file: pathlib.Path,
    attached_atom_idx: int,
    dummy_atom_idx: int,
) -> None:
    cube_data = flour.read_cube(cube_file)

    cube_data_positions_idx = _convert_euclidean_positions_to_indices(
        cube_data
    )

    streusel_surface = _get_surface(cube_data.grid.voxels)

    streusel_surface_idx = np.argwhere(streusel_surface)
    streusel_surface_point_cloud = pv.PolyData(streusel_surface_idx)

    p = pv.Plotter()
    p.add_mesh(streusel_surface_point_cloud)
    p.add_mesh(
        pv.PolyData(cube_data_positions_idx[dummy_atom_idx]),
        point_size=20,
        color="#69FAAB",
    )
    p.add_mesh(
        pv.PolyData(cube_data_positions_idx[attached_atom_idx]),
        point_size=20,
        color="#69FAAB",
    )
    p.show()


def _get_surface(voxels: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    cv = 0.0004
    is_vacuum = voxels < cv
    is_non_vacuum = np.logical_not(is_vacuum)

    weights = np.zeros((3, 3, 3))
    weights[1, 1, 2] = 1
    weights[1, 1, 0] = 1
    weights[:, :, 1] = [
        [0, 1, 0],
        [1, 0, 1],
        [0, 1, 0],
    ]

    convolution_result = scipy.ndimage.convolve(
        input=is_vacuum.view(np.int8),
        weights=weights,
        mode="wrap",
    )
    np.multiply(
        convolution_result,
        is_non_vacuum.view(np.int8),
        out=convolution_result,
    )
    return (convolution_result > 0).view(np.int8)


def _convert_euclidean_positions_to_indices(
    cube_data,  # : flour.CubeData,
) -> npt.NDArray:
    # translate origin and positions to 0,0,0
    translated_position_matrix = cube_data.positions - cube_data.grid.origin

    lattice_lengths = (
        cube_data.grid.voxel_size.sum(axis=0) * cube_data.grid.voxels.shape
    )

    percent_along_lattice_vector = np.divide(
        translated_position_matrix, lattice_lengths
    )

    return (percent_along_lattice_vector * cube_data.grid.voxels.shape).astype(
        int
    )


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Calculate and SMORES and sterimol parameters.",
    )

    parser.add_argument(
        "-i",
        "--input_file",
        help=(
            'A csv file with columns: "name", "core", "substituent", '
            '"smiles", "xyz_file", "dummy_index" and "attached_index".'
        ),
        type=pathlib.Path,
        default=pathlib.Path(
            "/home/bogosort/data/SMORES/implement-atomlite/validation/catalytic_reactions/one_substituent"
        )
        .joinpath("3_output")
        .glob("*/xyz_files.csv"),
        nargs="+",
    )
    parser.add_argument(
        "-o",
        "--output_directory",
        type=pathlib.Path,
        default=pathlib.Path.cwd() / "4_output",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
