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
