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

``generate_validation_systems.py`` generates the validation
system directories and initial structures.


## catalytic reactions
1. NHK allylations of benzaldehyde and acetophenone
2. desymmetrization of bisphenols
3. propargylation of acetophenone
