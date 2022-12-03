Validation with common carbon substituents
==========================================

In here you will find all the code required to compare
``SMORES`` parameters to sterimol parameters for
systems with common carbon substituents. You will find a script
for each step of the validation process, each helpfully beginning with
the step number. Each step can be run as-is, without any
configuration options, for example

.. code-block:: bash

  ./1_generate_systems
  ./2_calculate_electrostatic_potentials
  ./3_calculate_steric_parameters
  ./4_plot_steric_parameters

Each step produces an output folder which holds the input for
the following step. You can re-run individual steps,
provided the output folder for the previous step exists,
or the necessary input files are provided though options to the
script.

The difference between this validation and
common_carbon_substituents_fast__ is that the XYZ coordinates of the
tested structures are generated using ``Psi4`` rather than
``ETKDG``, making this validation significantly slower.

__ ../common_carbon_substituents_fast

To get more details about what each step does, or look at its
configuration options, use the ``-h`` parameter, for example

.. code-block:: bash

  ./1_generate_systems -h

The output of the final step will produce a bunch of nice graphs
for you to explore, which compare the calculated steric parameters
between different types of atomic radii and substituents.

Since each step produces a folder holding its output, you will end up
with the folders

* ``1_output``
* ``2_output``
* ``3_output``
* ``4_output``

if you run all your steps successfully.

Requirements
............

The ``conda`` environment file which has all the dependencies needed
to run the validation steps can be found here__.

__ ../smores.yml
