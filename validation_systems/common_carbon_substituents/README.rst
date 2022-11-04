Validation with common carbon substituents
==========================================

In here you will find all the code required to compare
``SMORES`` parameters to sterimol parameters for
systems with common carbon substituents. You will find a script
for each step of the validation process, each helpfully beginning with
the step number.

The difference between this validation and
common_carbon_substituents_fast__ is that the XYZ coordiantes of the
tested structures are generated using ``psi4`` rather than
``ETKDGv2``, making this validation significantly slower.

__ ../common_carbon_substituents_fast

To get more details about what each step does, use the ``-h``
parameter, for example

.. code-block:: bash

  ./1_generate_systems -h

Note that this may take some time, but will print output of the
following form::

  NAME: C_H
    SMORES L: ...   STERIMOL L: ...   Diff: ...
    SMORES B1: ...  STERIMOL B1: ...  Diff: ...
    SMORES B5: ...  STERIMOL B5: ...  Diff: ...
    RESULT: SUCCESS

  Total successful validations: X
  Total failed validations: Y

Each step will also produce a folder holding its output, you will end up
with the folders

* ``1_output``
* ``2_output``
* ``3_output``

Requirements
............

The ``conda`` envionrment file which has all the dependencies needed
to run the validation steps can be found here__.

__ ../../smores.yml
