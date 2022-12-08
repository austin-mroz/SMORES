#!python
import argparse
import csv
import pathlib
import typing
from dataclasses import dataclass

import ase.io.cube
import morfeus
import numpy as np

import smores
from smores._internal.write_cube import write_cube

_OUTPUT_CSV_COLUMNS = (
    "name",
    "core",
    "substituent",
    "smiles",
    "xyz_file",
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

                smores_molecule = smores.Molecule.from_xyz_file(row.xyz_file)
                smores_params = smores_molecule.get_steric_parameters(
                    dummy_index=row.dummy_index,
                    attached_index=row.attached_index,
                )
                writer.writerow(
                    {
                        "name": row.name,
                        "core": row.core,
                        "substituent": row.substituent,
                        "smiles": row.smiles,
                        "xyz_file": row.xyz_file,
                        "dummy_index": row.dummy_index,
                        "attached_index": row.attached_index,
                        "radii_type": "streusel",
                        "L": smores_params.L,
                        "B1": smores_params.B1,
                        "B5": smores_params.B5,
                    }
                )
                smores_esp_molecule = smores.EspMolecule.from_cube_file(
                    row.xyz_file.parent / "ESP.cube"
                )
                esp_smores_params = smores_esp_molecule.get_steric_parameters(
                    dummy_index=row.dummy_index,
                    attached_index=row.attached_index,
                )
                writer.writerow(
                    {
                        "name": row.name,
                        "core": row.core,
                        "substituent": row.substituent,
                        "smiles": row.smiles,
                        "xyz_file": row.xyz_file,
                        "dummy_index": row.dummy_index,
                        "attached_index": row.attached_index,
                        "radii_type": "streusel_cube",
                        "L": esp_smores_params.L,
                        "B1": esp_smores_params.B1,
                        "B5": esp_smores_params.B5,
                    }
                )
                data, atoms = ase.io.cube.read_cube_data(
                    str(row.xyz_file.parent / "ESP.cube")
                )
                vectors = np.linalg.norm(atoms.cell, axis=1)
                with open(row.xyz_file.parent / "ESP.cube") as esp_file:
                    lines = esp_file.readlines()
                write_cube(
                    path=args.output_directory
                    / f"{row.name}_{row.core}_{row.substituent}_esp.cube",
                    voxels=data,
                    voxel_dimensions=vectors / data.shape,
                    voxel_origin=np.array(
                        [float(i) for i in lines[2].split()[1:]]
                    ),
                    elements=atoms.get_atomic_numbers(),
                    positions=atoms.positions / ase.units.Bohr,
                )
                write_cube(
                    path=args.output_directory
                    / f"{row.name}_{row.core}_{row.substituent}.cube",
                    voxels=smores_esp_molecule._electric_field_surface.voxels,
                    voxel_dimensions=(
                        smores_esp_molecule._electric_field_surface.voxel_size
                    ),
                    voxel_origin=(
                        smores_esp_molecule._electric_field_surface.voxel_origin  # noqa
                    ),
                    elements=smores_esp_molecule._atoms,
                    positions=smores_esp_molecule._positions,
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
                            "dummy_index": row.dummy_index,
                            "attached_index": row.attached_index,
                            "radii_type": radii_type,
                            "L": sterimol.L_value,
                            "B1": sterimol.B_1_value,
                            "B5": sterimol.B_5_value,
                        },
                    )


@dataclass(frozen=True, slots=True)
class CsvRow:
    name: str
    core: str
    substituent: str
    smiles: str
    xyz_file: pathlib.Path
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
                dummy_index=int(row["dummy_index"]),
                attached_index=int(row["attached_index"]),
            )


def _get_sterimol(molecule: CsvRow, radii_type: str) -> morfeus.Sterimol:
    elements, coordinates = morfeus.read_xyz(molecule.xyz_file)
    return morfeus.Sterimol(
        elements=elements,
        coordinates=coordinates,
        dummy_index=molecule.dummy_index + 1,
        attached_index=molecule.attached_index + 1,
        radii_type=radii_type,
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
        default=pathlib.Path.cwd()
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
