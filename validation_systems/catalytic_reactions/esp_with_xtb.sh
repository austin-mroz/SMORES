#!/bin/bash

run_xtb () {

    cd $2
    DIR="$(dirname "$1")"
    cd "$2$DIR"
    xyz="initial_structure.xyz"
    xtb $xyz --esp
    cd "$2"
}

export -f run_xtb
parentDIR="/home/he/work/validation_calculations_for_SMORES_Oct_test_xtb/"
cd $parentDIR
find . -name "initial_structure.xyz" -exec bash -c 'run_xtb "$@" "'${parentDIR}'"' bash {} \;
