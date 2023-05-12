#!/bin/bash --login

#SBATCH --job-name=initialize-cages
#SBATCH --cpus-per-task=5
#SBATCH --time=48:30:00

#== Please do not change this section ==#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#======================================#

export LD_LIBRARY_PATH=/home/amroz/miniconda3/lib/

conda activate /home/amroz/SMORES/gen-psi4-esp/.venv

#### python 1_generate_systems.py
python 2_calculate_electrostatic_potentials.py
python 3_calculate_steric_parameters.py
python 4_plot_steric_parameters.py

