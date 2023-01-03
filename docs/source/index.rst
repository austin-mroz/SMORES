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
  rdkit_from_smiles <_autosummary/smores.rdkit_from_smiles>
  combine <_autosummary/smores.combine>
  psi4.calculate_electrostatic_potential <_autosummary/smores.psi4.calculate_electrostatic_potential>
  xtb.optimize_geometry <_autosummary/smores.xtb.optimize_geometry>
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

Getting help
------------

We want to make sure you're able to use :mod:`smores` to the
fullest, and we recognize that not everyone who wishes to use
:mod:`smores` is a confident programmer. As such, if you get
stuck using our tool we encourage you to get in touch with us on
Discord (`invite link`_), or by asking us in the `Q&A`_
section. We're happy to help!

If on the other hand you find an issue or bug with
:mod:`smores` please let us know by making an issue_.

.. _`invite link`: https://discord.gg/zbCUzuxe2B
.. _`Q&A`: https://github.com/austin-mroz/SMORES/discussions/categories/q-a
.. _issue: https://github.com/austin-mroz/SMORES/issues


Quickstart
----------

Getting started with :mod:`smores` is really simple!

You start with

.. doctest:: quickstart

  >>> import smores

and then you load a molecule and calculate the steric parameters


.. doctest:: quickstart

  >>> molecule = smores.Molecule.from_smiles("CC", dummy_index=0, attached_index=1)
  >>> molecule.get_steric_parameters()
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
library like sklearn_! Here we calculate the steric parameters for a
bunch of molecules with varying chain lengths and use them to predict
their UFF energy.

.. _sklearn: https://scikit-learn.org/stable/

.. testcode:: ml-workflow

  import rdkit.Chem.AllChem as rdkit
  import smores
  from sklearn.linear_model import LinearRegression

  def uff_energy(molecule):
      rdkit.SanitizeMol(molecule)
      return rdkit.UFFGetMoleculeForceField(molecule).CalcEnergy()

  cores = [smores.rdkit_from_smiles("CBr")]
  chains = ["C" * chain_length for chain_length in range(1, 50)]
  substituents = [smores.rdkit_from_smiles("Br" + chain) for chain in chains]

  X = []
  y = []
  for combo in smores.combine(cores, substituents):
      molecule = smores.Molecule.from_combination(combo)
      X.append(list(molecule.get_steric_parameters()))
      y.append(uff_energy(combo.product))

  reg = LinearRegression()
  reg.fit(X, y)

.. doctest:: ml-workflow

  >>> reg.score(X, y)
  0.8067966147933283

We hope that's a useful jumping off point for some quick prototyping!
We expect there are much more useful and interesting properties you can
target.

Plays nice with :mod:`rdkit`
----------------------------

:mod:`smores` molecules can easily be created from RDKit_ molecules

.. testcode:: rdkit-workflow

  import smores
  import rdkit.Chem.AllChem as rdkit

  rdkit_molecule = rdkit.AddHs(rdkit.MolFromSmiles("CBr"))
  rdkit.EmbedMolecule(rdkit_molecule)  # Generate a 3-D structure.
  smores_molecule = smores.Molecule.from_rdkit(rdkit_molecule, dummy_index=0, attached_index=1)


and we provide a handy function for creating rdkit
molecules from SMILES

.. testcode:: rdkit-workflow

   rdkit_molecule = smores.rdkit_from_smiles("CC")

which takes care of adding hydrogen atoms and generating a 3-D structure
for you, unlike RDKit_'s own `rdkit.MolFromSmiles`_ function.

.. _RDKit: https://www.rdkit.org/docs/index.html
.. _`rdkit.MolFromSmiles`: https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html#rdkit.Chem.rdmolfiles.MolFromSmiles


.. seealso::

   * :meth:`.Molecule.from_rdkit`: For addtional configuration options.
   * :meth:`.EspMolecule.from_rdkit`: For addtional configuration options.
   * :func:`.rdkit_from_smiles` For additional documentation.

Quick comparison of substituents
--------------------------------

A very common workflow is to try different substituents on a molecule
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
  for combo in smores.combine(cores, substituents):
      molecule = smores.Molecule.from_combination(combo)
      params = molecule.get_steric_parameters()
      print(
          f"Combination of {rdkit.MolToSmiles(rdkit.RemoveHs(combo.core))} and "
          f"{rdkit.MolToSmiles(rdkit.RemoveHs(combo.substituent))} "
          f"has SMORES parameters of {params}."
      )


.. testoutput:: substituent-comparison

  Combination of Brc1ccccc1 and CCCBr has SMORES parameters of StericParameters(L=5.6397512133212935, B1=1.7820154803719914, B5=3.4938688496917782).
  Combination of Brc1ccccc1 and CC(C)(C)CBr has SMORES parameters of StericParameters(L=5.668756954899209, B1=1.747631456476209, B5=4.532148595320116).

.. seealso::

   * :func:`.combine`: For additional examples and configuration options.
   * :func:`.rdkit_from_smiles`: For additional configuration options.


Using electrostatic potentials
------------------------------

:mod:`smores` can also calculate the steric parameters using electrostatic
potentials defined on a voxel grid


.. testcode:: using-electrostatic-potentials
  :hide:

  import os
  import tempfile
  tmp_dir = tempfile.TemporaryDirectory()
  os.chdir(tmp_dir.name)

.. doctest:: using-electrostatic-potentials

  >>> import smores
  >>> molecule = smores.EspMolecule.from_cube_file("HBr.cube", dummy_index=0, attached_index=1)
  >>> molecule.get_steric_parameters()
  StericParameters(L=3.57164113574581, B1=1.9730970556668774, B5=2.320611610648539)

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
  :hide:

  import os
  import tempfile
  tmp_dir = tempfile.TemporaryDirectory()
  os.chdir(tmp_dir.name)

.. testcode:: calculate-electrostatic-potential

  import smores
  import smores.psi4

  cube_path = smores.psi4.calculate_electrostatic_potential(
      molecule=smores.rdkit_from_smiles("Br"),
      output_directory="outdir",
      grid_origin=(-3., -3., -3.),
      grid_length=10.,
      num_voxels_per_dimension=20,
  )
  esp_molecule = smores.EspMolecule.from_cube_file(cube_path, dummy_index=0, attached_index=1)

.. doctest:: calculate-electrostatic-potential

  >>> esp_molecule.get_steric_parameters()
  StericParameters(L=3.57164113574581, B1=1.9730970556668774, B5=2.320611610648539)

.. seealso::

   * :func:`.psi4.calculate_electrostatic_potential`: For configuration options.

Optimizing molecules with xtb
-----------------------------

The values of the calculated steric parameters depend on the
geometry of the molecule. There are many different ways to
determine the geometry of a molecule. It may well be that the
structure returned by :func:`.rdkit_from_smiles` is suitable
for your purposes.

However, if that's not the case, you may want to optimize the
geometry yourself. To make your experience using
:mod:`smores` a little smoother, we provide a little helper
function in case you want to optimize your structures
using xtb_

.. _xtb: https://github.com/grimme-lab/xtb

.. testcode:: xtb-optimize-the-geometry-of-a-molecule
  :hide:

  import os
  import tempfile
  tmp_dir = tempfile.TemporaryDirectory()
  os.chdir(tmp_dir.name)

.. doctest:: xtb-optimize-the-geometry-of-a-molecule

  >>> import smores
  >>> molecule = smores.rdkit_from_smiles("CBr")
  >>> optimized = smores.xtb.optimize_geometry(molecule, "xtb_output")
  >>> smores_molecule = smores.Molecule.from_rdkit(optimized, dummy_index=0, attached_index=1)
  >>> smores_molecule.get_steric_parameters()
  StericParameters(L=3.57164113574581, B1=1.9730970556668774, B5=2.320611610648539)

Whether this is necessary or desirable, we leave to your judgement, but it's
here if you need it.

.. note::

 This example assumes that you have xtb available in your PATH.
 If that's not the case, you can set the location of the
 xtb binary using the `xtb_path` parameter. For more details,
 see :func:`.xtb.optimize_geometry`.

.. seealso::

   * :func:`.xtb.optimize_geometry`: For configuration options.
