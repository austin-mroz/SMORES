# overview
While we are waiting on the catalyst ESP's to be generated with Psi4, we will
work on integrating SMORES with STREUSEL(surface assessment) & 
morfeus/DBstep(sterimol parameter assessment).

We should also benchmark xtb ESP with Psi4 ESP for sterimol generation, as well.
^ to this end, have generated the xtb ESP gen submission script, but waiting to
run until the Psi4 ESPs are calculated for all of the full catalyst validation
systems.

# to do
- [X] export streusel neutral atom radii as a dictionary
- [X] insert streusel neutral atom radii dictionary in `src/smores/utiltiies.py`
- [X] generate list of streusel radii for each element in molecule
- [X] pipe into morfeus
- [X] generate common carbon substituent sterimol parameters for all radii types
offered in morfeus and sterimol >
`common_carbon_morfeus_xyz_sterimol_parameters.csv`
- [ ] generate sterimol parameters from streusel surface using cube file
- [ ] generate catalyst sterimol parameters for all radii types offered in
morfeus and sterimol
- [ ] algorithm to determine the atom indices for the sterimol calculations

### notes on the algorithm to determine the atom indices
In order to do this automatically, we should probably refactor the smores
organization because when the molecule is built, we should preserve the 'core'
and 'substituent' piece. This will allow us to save the indices that we will
need to calculate the sterimol parameters.
