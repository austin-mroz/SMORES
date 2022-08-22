# STREUSEL-based steric parameters
Essentially an article expanding the utility of the electric field size metric by allowing the calculation of additional steric
parameters at the electric field surface
To accomlish this we will develop a multi-step workflow, as follows:
1. 

# Outline
Here, we describe the implementation of varying steric metrics evaluated at the STREUSEL surface. This work merely presents an
extension to the STREUSEL package that requires minimal additional calculations to validate. We just have to show that we
perform similarly to accepted steric metrics.

# Relevant literature
| title | DOI |
| ----- | --- |
| Multidimensional steric parameters in the analysis of asymmetric catalytic reactions | 10.1038/NCHEM.1297 |
| | |
| Selection of low-dimensional 3D geometric descriptors for accurate enantioselectivity prediction | 10.1021/acscatal.2c00976 |
| | |
| Conformational effects on physical-organic descriptors: the case of Sterimol steric parameters | 10.1021/acscatal.8b04043 |

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

