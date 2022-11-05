.. smores documentation master file, created by
   sphinx-quickstart on Fri Nov  4 23:08:40 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to SMORES!
==================

.. toctree::
   :maxdepth: 2
   :caption: Contents:


GiHub: https://github.com/austin-mroz/SMORES


:mod:`smores` is a Python library which


It can lead to much better result than traditional
steric metrics, because ... See the paper ...



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

or directy from atomic coordinates

.. testcode:: quickstart

  molecule = smores.Molecule(
      atoms=["He"],
      positions=[[0., 0., 0.]],
  )

You can get the SMORES steric parameters by running

.. testcode:: quickstart

  params = molecule.get_steric_parameters()
  print(
      params.L,
      params.B1,
      params.B5,
  )

Which will calculate the parameters using the STREUSEL__ radii
of the atoms.

__ https://streusel.readthedocs.io


.. seealso::

  * :class:`.Molecule`
  * :meth:`.Molecule.get_steric_parameters`

More accurate parameters through electrostatic potentials
.........................................................

:mod:`smores` can also calculate the steric parameters using electrostatic
potentials defined on a voxel grid

.. testcode:: quickstart

  molecule = smores.EspMolecule.from_cube_file("my_molecule.cube")


.. testcode:: quickstart

  params = molecule.get_steric_parameters()
  print(
      params.L,
      params.B1,
      params.B5,
  )

.. seealso::

  * :class:`.EspMolecule`
  * :meth:`.EspMolecule.get_steric_parameters`

Installation
------------


.. code-block:: bash

  pip install smores


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
