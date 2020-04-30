:orphan:

.. raw:: html

  <div class="container-fluid">
    <div class="row">
      <div class="col-md-10">
        <div style="text-align: justify; text-justify: inter-word;">

Generate a Molecule
===================

This section describes the basics of creating molecule objects and their structures.


Initialize object
-----------------

An empty Molecule object can be easily created by passing no input. Remember to set a name to the molecule, which is also defined in the molecule config file.

.. code-block:: python

  from porems.molecule import Molecule

  mol = Molecule("benzene", "BEN")


Coordinate placement
--------------------

It is possible to create the structures by passing a coordinate.

.. code-block:: python

  mol.add("C", [0, 0, 0])


Position dependent placement
----------------------------

If the molecule has at least one atom, other atoms can be added depending on that position by passing the reference atom index a distance and an angle. In this case the angle depends on the z plane.

.. code-block:: python

  mol.add("C", 0, r=0.1375, theta=60)

.. figure::  /pics/mol/mol_1.png
  :align: center
  :width: 20%
  :name: fig2


Final molecule
--------------

Complete the ring structure of the benzene molecule.

.. code-block:: python

  mol.add("C", 1, r=0.1375, theta=120)
  mol.add("C", 2, r=0.1375, theta=180)
  mol.add("C", 3, r=0.1375, theta=240)
  mol.add("C", 4, r=0.1375, theta=300)

.. figure::  /pics/mol/mol_2.png
  :align: center
  :width: 30%
  :name: fig3


Save Structure
--------------

Using the storage class, the generated structure can be exported into various formats.

.. code-block:: python

  from porems.store import Store

  Store(mol).gro()
  Store(mol).pdb()
  Store(mol).xyz()


.. raw:: html

        </div>
      </div>
    </div>
  </div>
