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
  Modules <modules>

GiHub: https://github.com/austin-mroz/SMORES


:mod:`smores` is a Python library which


It can lead to much better result than traditional
steric metrics, because ... See the paper ...

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

and then you load a molecule, either from one of the supported
file formats

.. testcode:: quickstart

  molecule = smores.Molecule.from_xyz_file("my_molecule.xyz")

or from SMILES

.. testcode:: quickstart

  molecule = smores.Molecule.from_smiles("CBr")

or directy from atomic coordinates

.. testcode:: quickstart

  molecule = smores.Molecule(
      atoms=["He"],
      positions=[[0., 0., 0.]],
  )

You can get the SMORES steric parameters by running

.. testcode:: quickstart

  params = molecule.get_steric_parameters()
  print(params.L, params.B1, params.B5)

Which will calculate the parameters using the STREUSEL__ radii
of the atoms.

__ https://streusel.readthedocs.io


.. seealso::

  * :class:`.Molecule`
  * :meth:`.Molecule.get_steric_parameters`

Loading with atomic numbers
---------------------------

Sometimes you have atomic numbers but not
atomic elements, fortunately we have you covered

.. testcode:: quickstart

  molecule = smores.Molecule(
      atoms=smores.atomic_numbers_to_elements([1, 35]),
      positions=[[0., 0., 0.], [1.47, 0., 0.]],
  )

Using electrostatic potentials
------------------------------

:mod:`smores` can also calculate the steric parameters using electrostatic
potentials defined on a voxel grid

.. testcode:: quickstart

  molecule = smores.EspMolecule.from_cube_file("my_molecule.cube")
  params = molecule.get_steric_parameters()
  print(params.L, params.B1, params.B5)

.. tip::

  There are more ways to initialize an :class:`.EspMolecule`! See the
  documentation of the class for details.

.. seealso::

  * :class:`.EspMolecule`
  * :meth:`.EspMolecule.get_steric_parameters`






Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
