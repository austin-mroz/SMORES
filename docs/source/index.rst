.. smores documentation master file, created by
   sphinx-quickstart on Fri Nov  4 23:08:40 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to SMORES!
==================

.. toctree::
  :hidden:
  :maxdepth: 2
  :caption: API

  Molecule <_autosummary/smores.Molecule>
  EspMolecule <_autosummary/smores.EspMolecule>
  combine <_autosummary/smores.combine>
  Psi4 <_autosummary/smores.psi4>
  Modules <modules>

GitHub: https://github.com/austin-mroz/SMORES


.. attention::

  This introduction is under construction but will be
  filled in by Austin + Lukas.

  Things to cover:

  * What does smores do?
  * Why would the user want to do this?
  * What are steric metrics, why are they useful?
  * What can you do with steric metrics?
  * Why is using eletric field based stuff cool?
  * Link to STREUSEL
  * A cool graphic never hurt anyone
  * Talk about morfeus


Installation
------------


.. code-block:: bash

  pip install smores


Quickstart
----------

Getting started with :mod:`smores` is really simple!

You start with

.. testcode:: quickstart

  import smores

and then you load a molecule & calculate the steric parameters


.. testcode:: quickstart

  molecule = smores.Molecule.from_smiles("CC")
  molecule.get_steric_parameters(dummy_index=0, attached_index=1)

.. testoutput::

  StericParameters(L=3.57164113574581, B1=1.9730970556668774, B5=2.320611610648539)


Which will calculate the parameters using the STREUSEL__ radii
of the atoms.

__ https://streusel.readthedocs.io

.. tip::

  The radii which are used can be modified by the user! See
  the documentation of :meth:`.Molecule.get_steric_parameters`
  for more details.

.. seealso::

  * :class:`.Molecule`: For additional documentation and examples.
  * :meth:`.Molecule.from_xyz_file`: For loading molecules from ``.xyz`` files.
    Other file types are supported too!
  * :meth:`.Molecule.get_steric_parameters`: For configuration options.

Integration with machine learning workflows
-------------------------------------------

It doesn't take a lot of code to get :mod:`smores` working with a great
library like sklearn_

.. _sklearn: https://scikit-learn.org/stable/

.. testcode:: ml-workflow

  import smores
  import sklearn
  from glob import glob

  molecules = [smores.Molecule.from_xyz_file(path) for path in glob("*.xyz")]

  # An N x 3 array, where each row holds L, B1 and B2 of a given molecule.
  params = np.array([list(mol.get_steric_parameters(0, 1)) for mol in molecules])

  classifier = sklearn.tree.DecisionTreeClassifier()
  classifier.fit(params, target_property)
  classifier.predict(smores.Molecule.from_xyz_file("test.xyz"))


We hope that's a useful jumping off point for some quick prototyping!

Plays nice with :mod:`rdkit`
----------------------------

:mod:`smores` molecules can easily be created from RDKit_ molecules

.. testcode::

  import smores
  import rdkit.Chem as rdkit

  rdkit_molecule = rdkit.AddHs(rdkit.MolFromSmiles("CBr"))
  rdkit.EmbedMolecule(rdkit_molecule)  # Generate a 3-D structure.

  smores_molecule = smores.Molecule.from_rdkit(rdkit_molecule)


and we provide a handy function for creating rdkit
molecules from SMILES

.. testcode::

   rdkit_molecule = smores.rdkit_from_smiles("CC")

which takes care of adding hydrogen atoms and generating a 3-D structure
for you, unlike RDKit_'s own `rdkit.MolFromSmiles`_ function.

.. _RDKit: https://www.rdkit.org/docs/index.html
.. _`rdkit.MolFromSmiles`: https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html#rdkit.Chem.rdmolfiles.MolFromSmiles


.. seealso::

   * :meth:`.Molecule.from_rdkit`: For addtional configuration options.
   * :meth:`.EspMolecule.from_rdkit`: For addtional configuration options.

Quick comparison of substituents
--------------------------------

A very common workflow is to try different subsituents on a molecule
and compare their steric parameters, so we wrote some code that lets
you do this quick

.. testcode:: substituent-comparison

  import smores
  import rdkit.Chem as rdkit

  cores = [
      smores.rdkit_from_smiles("c1ccccc1Br"),
  ]
  substituents = [
    smores.rdkit_from_smiles("BrCCC"),
    smores.rdkit_from_smiles("BrCC(C)(C)C"),
  ]
  for combo in smores.combine(cores, subsituents):
      molecule = smores.Molecule.from_rdkit(combo.product)
      params = molecule.get_steric_parameters(
          dummy_index=combo.dummy_index,
          attached_index=combo.attached_index,
      )
      print(
          f"Combination of {rdkit.MolToSmiles(combo.core)} and "
          f"{rdkit.MolToSmiles(combo.substituent)} "
          f"has SMORES parameters of {params}."
      )


Note that this code allows you to easily identify what the dummy
and attachment atoms are, which can be a bit of a burden otherwise.


.. seealso::

   * :func:`.combine`: For additional examples and configuration options.
   * :func:`.rdkit_from_smiles`: For additional configuration options.


Using electrostatic potentials
------------------------------

:mod:`smores` can also calculate the steric parameters using electrostatic
potentials defined on a voxel grid

.. testcode:: quickstart

  molecule = smores.EspMolecule.from_cube_file("my_molecule.cube")
  params = molecule.get_steric_parameters(dummy_index=0, attached_index=1)
  print(params.L, params.B1, params.B5)


.. seealso::

  * :class:`.EspMolecule`: For additional documentation and examples.
  * :meth:`.EspMolecule.get_steric_parameters`: For configuration options.
  * :mod:`smores.psi4`: For using Psi4 to make ``.cube`` files.

Calculating electrostatic potentials
------------------------------------

.. attention::

  Psi4_ is a big dependency and we therefore do not install it automatically.
  If you wish to use the functions in :mod:`smores.psi4` you will have to install
  Psi4_ yourself. We find

  .. code-block:: bash

    conda install -c psi4 psi4

  tends to work for us, but we recommend you check out Psi4_'s official
  documentation for up-to-date information.

.. _psi4: https://psicode.org/


.. testcode:: calculate-electrostatic-potential

  import smores
  import smores.psi4

  smores.psi4.calculate_electrostatic_potential(
      molecule=smores.rdkit_from_smiles("Br"),
      output_directory="outdir",
      grid_origin=(-3., -3., -3.),
      grid_length=10.,
      num_voxels_per_dimension=50,
  )

  esp_molecule = smores.EspMolecule.from_cube_file("outdir/ESP.cube")
  params = esp_molecule.get_steric_parameters(dummy_index=0, attached_index=1)
  print(params.L, params.B1, params.B5)

.. seealso::

   * :func:`.calculate_electrostatic_potential`: For configuration options.
