"""
This module provides utility functions for using Psi4.

.. attention::

    Psi4_ is a big dependency and we therefore do not install it automatically.
    If you wish to use the functions in this module you will have to install
    Psi4_ yourself. We find

    .. code-block:: bash

        conda install -c psi4 psi4

    tends to work for us, but we recommend you check out Psi4_'s official
    documentation for up-to-date information.

.. _psi4: https://psicode.org/

Examples:

    Use the electrostatic potential, to calculate the SMORES
    parameters

    .. testcode:: calculate-electrostatic-potential

        import smores
        import smores.psi4

        molecule = smores.Molecule.from_smiles("Br")
        esp = smores.psi4.calculate_electrostatic_potential(molecule, "outdir")
        esp_molecule = smores.EspMolecule(
            atoms=molecule.atoms,
            positions=molecule.positions,
            electrostatic_potential=esp,
        )
        params = esp_molecule.get_steric_parameters(dummy_index=0, \
attached_index=1)
        print(params.L, params.B1, params.B5)


"""


from smores._internal.psi4 import calculate_electrostatic_potential

__all__ = [
    "calculate_electrostatic_potential",
]
