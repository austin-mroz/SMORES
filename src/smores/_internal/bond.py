from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class Bond:
    atom1: int
    atom2: int
    order: int
