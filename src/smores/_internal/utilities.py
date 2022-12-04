import contextlib
import os
import pathlib
import typing


@contextlib.contextmanager
def current_working_directory(directory: pathlib.Path) -> typing.Any:
    original_directory = os.getcwd()  # aka OGD
    try:
        os.chdir(directory)
        yield
    finally:
        os.chdir(original_directory)
