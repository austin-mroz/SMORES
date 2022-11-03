import argparse
import csv
import itertools
import pathlib
import typing
from dataclasses import dataclass


def main() -> None:
    args = _get_command_line_arguments()
    args.output_file.parent.mkdir(parents=True, exist_ok=True)

    substituents = (
        Substituent(name="H", smiles="[DUMMY]"),
        Substituent(name="Me", smiles="[DUMMY]C"),
        Substituent(name="Et", smiles="[DUMMY]CC"),
        Substituent(name="Ph", smiles="[DUMMY]c1ccccc1"),
        Substituent(name="Bn", smiles="[DUMMY]Cc1ccccc1"),
        Substituent(name="CH2-iPr", smiles="[DUMMY]C(C)CC"),
        Substituent(name="CH2-tBu", smiles="[DUMMY]C(C)(C)CC"),
        Substituent(name="i-Pr", smiles="[DUMMY](C)CC"),
        Substituent(name="CHPr2", smiles="CCCC([DUMMY])CCC"),
        Substituent(name="Cy", smiles="N#C[DUMMY]"),
        Substituent(name="CH(i-Pr)2", smiles="CC(C)C([DUMMY])C(C)C"),
        Substituent(name="CHEt2", smiles="CCC([DUMMY])CC"),
        Substituent(name="CEt3", smiles="CCC([DUMMY])(CC)CC"),
    )
    r_groups = ("C", "c1ccccc1")

    with open(args.output_file, "w") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=["name", "smiles"])
        writer.writeheader()
        writer.writerows(_get_csv_rows(substituents, r_groups))


@dataclass(frozen=True, slots=True)
class Substituent:
    name: str
    smiles: str


def _get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate the SMILES of systems which are to be validated.",
    )
    parser.add_argument(
        "-o",
        "--output_file",
        default=pathlib.Path(__file__).parent / "1_output" / "smiles.csv",
        type=pathlib.Path,
        help="The csv file into which the generated SMILES are written.",
    )
    return parser.parse_args()


class CsvRow(typing.TypedDict):
    name: str
    smiles: str


def _get_csv_rows(
    substituents: typing.Iterable[Substituent],
    r_groups: typing.Iterable[str],
) -> typing.Iterator[CsvRow]:

    for substituent, r_group in itertools.product(substituents, r_groups):
        yield CsvRow(
            name=f"{r_group}_{substituent.name}",
            smiles=substituent.smiles.replace("[DUMMY]", r_group),
        )


if __name__ == "__main__":
    main()
