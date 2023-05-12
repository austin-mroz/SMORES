#!/bin/bash --login

#SBATCH --job-name=initialize-cages
#SBATCH --cpus-per-task=12
#SBATCH --time=24:30:00

#== Please do not change this section ==#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#======================================#

export LD_LIBRARY_PATH=/home/amroz/miniconda3/lib/

conda activate /home/amroz/SMORES/gen-psi4-esp/.venv

#### python 1_generate_systems.py
python 2_calculate_electrostatic_potentials.py
python 3_calculate_steric_parameters.py -o /home/amroz/SMORES/sterics-from-cubefile/validation/common_carbon_substituents/3_output_streusel_cube
python 4_plot_steric_parameters.py -i /home/amroz/SMORES/sterics-from-cubefile/validation/common_carbon_substituents/3_output_streusel_cube/steric_parameters.csv -o /home/amroz/SMORES/sterics-from-cubefile/validation/common_carbon_substituents/4_output_streusel_cube

