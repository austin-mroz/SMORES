Validation with catalytic systems
=================================

In here you will find all the code required to compare
``SMORES`` parameters to sterimol parameters for catalytic
systems described by Harper__ et al. You will find a script
for each step of the validation process, each helpfully beginning with
the step number. Each step can be run as-is, without any
configuration options, for example

__ https://www.nature.com/articles/nchem.1297

.. code-block:: bash

  ./1_generate_systems
  ./2_optimize_strucures
  ./3_calculate_electrostatic_potentials
  ./4_compare_smores_with_sterimol

Each step produces an output folder which holds the input for
the following step. You can re-run individual steps,
provided the output folder for the previous step exists,
or the necessary input files are provided though options to the
script.

To get more details about what each step does, or look at its
configuration options, use the ``-h`` parameter, for example

.. code-block:: bash

  ./1_generate_systems -h

The output of the final step will produce output of the following
form::

  NAME: C_H
    SMORES L: ...   STERIMOL L: ...   Diff: ...
    SMORES B1: ...  STERIMOL B1: ...  Diff: ...
    SMORES B5: ...  STERIMOL B5: ...  Diff: ...
    RESULT: SUCCESS

  Total successful validations: X
  Total failed validations: Y

Since each step produces a folder holding its output, you will end up
with the folders

* ``1_output``
* ``2_output``
* ``3_output``
* ``4_output``

if you run all your steps successfully.

Requirements
............

You must have xtb__ installed and in your ``PATH``.

__ https://xtb-docs.readthedocs.io/en/latest/contents.html
