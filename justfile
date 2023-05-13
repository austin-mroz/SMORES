# list available recipes
default:
  @just --list

# build the docs
docs:
  rm -rf docs/build docs/source/_autosummary
  make -C docs html
  echo Docs are in $PWD/docs/build/html/index.html

# apply formatters to the code
fix:
  black .
  isort .

# check the code for errors
check:
  #!/usr/bin/env bash

  error=0
  trap error=1 ERR

  echo
  ( set -x; black --check . )
  echo
  ( set -x; isort --check . )
  echo
  ( set -x; flake8 . )
  echo
  ( set -x; mypy --exclude validation . )
  echo
  ( set -x; mypy --explicit-package-bases validation )
  echo
  ( set -x; pytest )

  test $error = 0

# build the docker testing environment
build-testing-environment:
  pip-compile -o docker_testing_environment/requirements.txt --extra dev pyproject.toml
  docker image build -t smores-testing-environment:latest docker_testing_environment

# create python 3.10 environment
venv-310:
  conda create -p .py310 python=3.10

# create python 3.11 environment
venv-311:
  conda create -p .py311 python=3.11
