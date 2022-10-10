# validating-systems
Here, we outline the systems that we'll use to validate SMORES.
These are largely drawn from Sigman's Nat. Chem. paper

The systems are organized using the following structure:

```
validation-systems
|   validation_plan.md
|
|---common-carbon-substituents
|   |---H
|   |   *sub*
|   |---Me
|   |   *sub*
|   |---Bn
|   |   *sub*
|---catalytic-reactions
|   |---NHK
|   |---desymmetrization
|   |---propargylation
```
## common carbon substituents
We will generate molecular systems where the following groups
are substituents on a) H, b) CH3, and c) Bn
1. H
2. Me
3. Et
4. Ph
5. Bn
6. CH2-iPr
7. CH2-tBu
8. i-Pr
9. CHPr2
10. Cy
11. CH(i-Pr)2
12. CHEt2
13. CEt3

Was unable to identify the chemical structure of the final
substituent:
14. Ad

We will generate the substituent validation systems using a HT
approach. The substituents are defined in the
``carbon_substituents_smiles.csv``. Br is used as the
placeholder group for HT generation.

``generate_common_carbon_validation_systems.py`` generates the validation
system directories and initial structures.


## catalytic reactions
We will validate SMORES with the following catalytic reactions:
1. NHK allylations of benzaldehyde
2. NHK allylation of acetophenone
3. desymmetrization of bisphenols
4. propargylation of acetophenone

The catalysts for each of these are found in:
`catalytic_reaction_smiles.csv`
The subsitunt placement is denoted by `Br` for the cases where there is only one
substituent (1,2,3). In cases where 2 substiutents are used (1,4), `Lu` is used
as the second placement atom for enumeration purposes.

`generate_catalyst_validation_systems.py` generates the validation system
dictectories, initial structures, as well as optimizing the catalysts.

# progress
Common carbon substituent systems have been generated.

Catalytic reaction validation systems are partly completed --
Presently, the code can only handle catalysts with one substituent group
placement (assumed to be Br). However, both the NHK propargylation and the NHK
allylation with the oxazoline-proline-library require 2 substituent groups
(assumed to be Br and Lu).

`generate_catalyst_validation_systems.py` needs to be modified to accomodate
these cases.

should also add an option to optimize with xtb instead of psi4 -- xtb would be
much faster

