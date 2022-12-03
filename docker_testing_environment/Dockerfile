FROM python:3.10.8-slim

# pyproject.toml needs to be copied into the build context
# from the root of the repository.
COPY pyproject.toml /deps/pyproject.toml

RUN python -m venv /deps/venv && . /deps/venv/bin/activate && pip install '/deps/[dev]'
