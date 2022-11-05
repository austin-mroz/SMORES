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
  Psi4 <_autosummary/smores.psi4>
  Modules <modules>

GiHub: https://github.com/austin-mroz/SMORES


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

  molecule = smores.Molecule.from_xyz_file("my_molecule.xyz")
  params = molecule.get_steric_parameters(dummy_index=0, attached_index=1)
  print(params.L, params.B1, params.B5)

Which will calculate the parameters using the STREUSEL__ radii
of the atoms.

__ https://streusel.readthedocs.io

.. tip::

  The radii which are used can be modified by the user! See
  the documentation of :meth:`.Molecule.get_steric_parameters`
  for more details.

.. seealso::

  * :class:`.Molecule`: For additional documentation and examples.
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
      molecule=smores.Molecule.from_smiles("Br"),
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
