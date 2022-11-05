import typing
from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class StericParameters:
    """
    The SMORES steric parameters.

    """

    #: The L parameter.
    L: float

    #: The B\ :sub:`1` parameter.
    B1: float

    #: The B\ :sub:`5` parameter.
    B5: float

    def __iter__(self) -> typing.Iterator[float]:
        yield self.L
        yield self.B1
        yield self.B5
