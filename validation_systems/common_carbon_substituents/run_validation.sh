#!/usr/bin/env bash

python 1_generate_systems.py
python 2_calculate_electrostatic_potentials.py
python 3_compare_smores_with_sterimol.py
