Validation with common carbon substituents
==========================================

In here you will find all the code required to compare
``SMORES`` parameters to sterimol parameters for
systems with common carbon substituents. You will find a script
for each step of the validation process, helpfully beginning with the
step number.

To get more details about what each step does, use the ``-h``
parameter, for example

.. code-block:: bash

  python 1_generate_systems.py -h

You can run all the steps in order with

.. code-block:: bash

  ./run_validation.sh

Note that this may take some time, but will print output of the
following form::

  SMILES: C
    SMORES L: ...   STERIMOL L: ...   Diff: ...
    SMORES B_1: ... STERIMOL B_1: ... Diff: ...
    SMORES B_5: ... STERIMOL B_5: ... Diff: ...
    RESULT: SUCCESS

  Total successful validations: X
  Total failed validations: Y

Each step will also produce a folder holding it output, you will end up
with the folders

* ``1_output``
* ``2_output``
* ``3_output``
