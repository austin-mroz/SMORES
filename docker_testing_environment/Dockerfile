FROM python:3.10.8-slim

# requirements.txt is generated with pip-compile.
COPY requirements.txt /deps/requirements.txt

RUN python -m venv /deps/venv && . /deps/venv/bin/activate && pip install -r /deps/requirements.txt
