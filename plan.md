# STREUSEL-based steric parameters
Essentially an article expanding the utility of the electric field size metric by allowing the calculation of additional steric
parameters at the electric field surface

## functionality
The major functionalities are outlined by the multi-step workflow:
1. efficient electric field generation (2 potential avenues)
	1. electric field from multipole expansion -- where multipoles are predicted using [6]`https://github.com/rinikerlab/EquivariantMultipoleGNN`
	2. electric field from electrostatic potential predicted using [5]`https://github.com/AstexUK/ESP_DNN`
        1. ESP_DNN outputs .pqr files > by extension, this means that [5] actually predicts atomic charges and radii, which are the final two columns in a .pqr file
    3. this will be validated against the results presented in the original STREUSEL paper and DeRose-JACS implementation of STREUSEL. Visualization can be similar to that presented in the GNN (EquivariantMultipoleGNN paper)
2. electric field-based sterimol (via DBStep) and weighted sterimol (via wSterimol, [3]) parameters will be generated
	1. weighted and unweighted sterimol parameters (L, B1, B2) will be validated/compared with those presented by DBStep and wSterimol [3]
3. electric field-based sterimol parameters will be used to parameterize the catalytic examples in [1]



# Outline
Here, we describe the implementation of varying steric metrics evaluated at the STREUSEL surface. This work merely presents an
extension to the STREUSEL package that requires minimal additional calculations to validate. We just have to show that we
perform similarly to accepted steric metrics.

# Relevant literature
| Ref. No. | title | DOI |
| -------- | ----- | --- |
| 1 | Multidimensional steric parameters in the analysis of asymmetric catalytic reactions | 10.1038/NCHEM.1297 |
| | | |
| 2 | Selection of low-dimensional 3D geometric descriptors for accurate enantioselectivity prediction | 10.1021/acscatal.2c00976 |
| | | |
| 3 | Conformational effects on physical-organic descriptors: the case of Sterimol steric parameters | 10.1021/acscatal.8b04043 |
| | | |
| 4 | Practical High-Quality Electrostatic Potential Surfaces for Drug Discovery Using a Graph-Convolutional Deep Neural Network | 10.1021/acs.jmedchem.9b01129 |
| | | |
| 5 | ESP-DNN: a graph-convolutional DNN for predicting electrostatic potential surfaces | 10.1021/acs.jmedchem.9b01129 |
| | | |
| 6 | Learning atomic multipoles: prediction of the electrostatic potential with equivariant GNNs | 10.1021/acs.jctc.1c01021 |
| | | |


# Steric metrics we will compare against

## Classical steric parameters
These parameters are introduced by <mark>10.1038/NCHEM.1297</mark>

| parameter | definition | reference | 
| --------- | ---------- | --------- |
| Winstein-Holness | A-values, conformational study of mono-substituted cyclohexane rings | |
| | | |
| Molar refractivity | defined by Lorentz-Lorenz equation, disregards molecular shape | |
| | | |
| Tolman cone angle | likely limited to phosphine ligands | | 
| | | |
| Taft parameter | | |
| | | |
| Charton | correlation between Taft experimental rates and min vdW radii of each symmetrical substituent | |
| | | |
| Hansch | validation of Charton and extrapolates Charton to unmeasured substituents | |
| | | |

# Notes
Perhaps we consider a STREUSEL extension to DBStep or wSterimol, which were developed by Robert Paton.

First step is to implement dbstep in streusel -- then we can more readily determine whether this will work
on a larger scale. We will start with looking at ethane, a simple example molecule that is proposed on the 
DBSTEP github.

To implement sterimol in streusel through dbstep, we will import streusel radii into dbstep (which already
has an implementation for bondi radii). Thus, it is a simple modification to add a tag calling the streusel
radii.
Then, we will examine the use of streusel surface. To do this the simplest way possible, we will take
advantage of the surface feature of dbstep, but instead of calling a normal cube file, we will rewrite a cube
file that denotes the streusel surface as 0.002 isoval, which is the magnitude that dbstep looks for in its
program.
Should we obtain realistic/comparable sterimol values, we will consider full implementation in the streusel
package and further examine the correlations observed in the case studies found in the 10.1021/acscatal.8b04043
paper cited above (wSterimol paper presented by Paton).

# Benchmark progress
We have generated a cube file for the ethane.xyz structure presented on the DBstep github page.
We will validate SMORES via several systems. These are located:
`/SMORES/validation-systems/*/`. Scripts that were used to generate the systems can also be found in this
directory.

# Project progress
We have generated the systems that will be used to validate SMORES.
The first 2 steps of the 3 step worklow are complete.
1. ESP maps may be generated using Psi4
2. STREUSEL has been integrated with SMORES -- this reuqires additional validation steps (for the optimized
binning method).
3. SMORES sterimol paramters must yet be integrated with MOR 


