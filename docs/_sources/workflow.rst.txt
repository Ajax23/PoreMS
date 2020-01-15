:orphan:

.. raw:: html

  <br>

.. raw:: html

  <div class="container-fluid">
    <div class="row">
      <div class="col-md-10">
        <div style="text-align: justify; text-justify: inter-word;">

Simulation Workflow
===================

In this workflow a simply pore simulation system will be created with TMS as
surface molecules.


Create surface molecules
------------------------

Trimethysilyl or for short TMS is a simple surface group that can be imported
from the porems package

.. code-block:: python

  from porems.essentials import TMS

Following code block can be used as a base structure to create further surface
group structures

.. code-block:: python

  from porems.molecule import Molecule

  tms = Molecule("tms", "TMS")
  tms.set_charge(0.96)
  compress = 30

  b = {"sio": 0.161, "sic": 0.190}
  a = {"ccc": 30.00, "cch": 109.47}

  tms.add("Si", [0, 0, 0])
  tms.add("O", 0, r=b["sio"])
  tms.add("Si", 1, bond=[1, 0], r=b["sio"], theta=180)

  for i in range(3):
      tms.add("C", 2, bond=[2, 1], r=b["sic"], theta=a["cch"]+compress, phi=60+120*i)



Create pore system
------------------

.. code-block:: python

  from porems.pore import Pore

  pore = Pore(size=[10, 10, 10], diam=6, drill="z", res=5.5, is_time=True)

  pore.siloxan(0.5)

  pore.attach(tms, [0, 1], [1, 2], 0, 3, inp="molar", is_rotate=False)
  pore.attach(tms, [0, 1], [1, 2], 1, 3, inp="molar", is_rotate=False)

  pore.finalize()


.. code-block:: python

  Store(pore).gro("pore.gro")
  Store(pore).obj("pore.obj")
  Store(pore).top("topol.top")
  Store(pore).grid("grid.itp")


Simulation folder structure
---------------------------

Provide zip file of basic structure (nofill)


Fixiating surface molecules and grid
------------------------------------

gromacs index call


Filling box
-----------

Example water density of tip3p GAFF


Basic analysis idea of density
------------------------------

Formulas and import data from pore object


.. raw:: html

        </div>
      </div>
    </div>
  </div>
