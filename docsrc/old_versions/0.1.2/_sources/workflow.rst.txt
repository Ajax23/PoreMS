:orphan:

.. raw:: html

  <div class="container-fluid">
    <div class="row">
      <div class="col-md-10">
        <div style="text-align: justify; text-justify: inter-word;">

Simulation Workflow
===================

In this workflow a simply pore simulation system will be created with TMS as
surface molecules. Additionally, the GROMACS simulation package will be utilized.


Create surface molecules
------------------------

Trimethylsilyl or for short TMS is a simple surface group that can be imported
from the porems package.

.. code-block:: python

  from porems.essentials import TMS

Assuming a new surface group structure is to be created, following code block
can be used as a base.

.. code-block:: python

  from porems.molecule import Molecule

  tms = Molecule("tms", "TMS")
  tms.set_charge(0.96)
  compress = 30

  b = {"sio": 0.155, "sic": 0.186, "ch": 0.109}
  a = {"ccc": 30.00, "cch": 109.47}

  # Create tail
  tms.add("Si", [0, 0, 0])
  tms.add("O", 0, r=b["sio"])
  tms.add("Si", 1, bond=[1, 0], r=b["sio"], theta=180)

  # Add carbons
  for i in range(3):
      tms.add("C", 2, bond=[2, 1], r=b["sic"], theta=a["cch"]+compress, phi=60+120*i)

  # Add hydrogens
  for i in range(3, 5+1):
      for j in range(3):
          tms.add("H", i, bond=[i, 2], r=b["ch"], theta=a["cch"]+compress, phi=60+120*j)

.. figure::  /pics/flow/tms.png
 :align: center
 :width: 30%
 :name: fig1

.. note::

  Parametrization must be carried out by the user. Topology generation should
  be performed for both a singular binding site and a geminal binding site.


Create pore system
------------------

Next step is to create a pore structure functionalized with the created TMS
surface group.

.. code-block:: python

  from porems.pore import Pore

  pore = Pore(size=[10, 10, 10], diam=6, drill="z", res=5.5, is_time=True)

  pore.siloxan(0.5)

  pore.attach(tms, [0, 1], [1, 2], 0, 3, inp="molar", is_rotate=False)
  pore.attach(tms, [0, 1], [1, 2], 1, 3, inp="molar", is_rotate=False)

  pore.finalize()

.. figure::  /pics/flow/pore.png
 :align: center
 :width: 50%
 :name: fig2

Once the generation is done, store the structure and preferably the object for
future analysis. Furthermore, a master topology with the number of residues and
a topology containing grid molecule parameters should be created using the
:func:`porems.store.Store.top` and :func:`porems.store.Store.grid` functions.

.. code-block:: python

  Store(pore).gro("pore.gro")
  Store(pore).obj("pore.obj")
  Store(pore).top("topol.top")
  Store(pore).grid("grid.itp")


Simulation folder structure
---------------------------

The simulations folder :download:`provided <data/test_sim.zip>` has following structure

* Top Folder

  * **_top** - Folder containing topologies

    * **topol.top** - Master topology

    * **grid.itp** - Grid molecule parameters

    * **tip3p.itp** - Topology for TIP3P water

    * **tms.itp** - Topology for TMS with singular binding site

    * **tmsg.itp** - Topology for TMS with geminal binding site

  * **_gro** - Folder containing structure files

    * **box.gro** - Simulation box

    * **spc216.gro** - Water structure to be filled in the simulation box

    * **pore.obj** - Pore object as a backup for future analysis

  * **_mdp** - Folder containing simulation parameter files

    * **min.mdp** - Energy minimization parameter file

    * **nvt.mdp** - NVT equilibration parameter file

    * **run.mdp** - Production parameter file

  * **min** - Folder for carrying out energy minimization

  * **nvt** - Folder for carrying out NVT equilibration

  * **run** - Folder for the production run

.. note::

  Topologies provided are from the General AMBER Force Field (GAFF).

  Furthermore, the excess charge which might arise from surface molecule
  parametrization can be distributed among the grid molecules.



Fixating surface molecules and grid
------------------------------------

The grid is fixated by removing specified atoms from the energy calculation of
GROMACS. This can be done by first defining an index group

.. code-block:: bash

  gmx make_ndx -f _gro/box.gro -o _gro/index.ndx

and choosing the specified atoms. Since ``make_ndx`` works iterativaly, first
the silicon atoms of the surface groups, silanol and TMS, are chosen for both
geminal and singular binding sites, and then the grid molecules. In the case of
the generated pore system, the call would be

.. code-block:: bash

  5 | 6 | 7 | 8 & a SI1 | 2 | 3 | 4

This index group is then specified in the mdp files under the freezed groups tag

.. code-block:: bash

  freezegrps = SL_SLG_TMS_TMSG_&_SI1_SI_OM_OX
  freezedim  = Y Y Y

.. note::

   To make sure all fixated atoms were added to the index group, a simple
   calculation should be performed before simulation.



Filling box
-----------

The pore system is simulated in the NVT ensample, since NPT in GROMACS displaces
the grid molecules in the simulation while adjusting the box-size to the pressure.
Nonetheless, the system needs to be simulated at a specified density. This is
done by iteratively filling the box with the solute molecules, here water, until
achieving the reference density as in an NPT simulation at the desired pressure.

.. note::

  If the GROMACS filling functions, like ``solvate`` or ``insert-molecules``
  are used with small molecules, it may happen that molecules are placed within
  the grid. Naturally these molecules must be removed from the grid before
  running the simulation.


Density analysis procedure
--------------------------

The density calculation inside and outside the pore is done by calculating
the number density :math:`\rho_n` and using the molar mass :math:`M` of the
molecule to determine the mass density :math:`\rho`.

The basic idea is counting the number of molecules :math:`N_i` in volume slices
:math:`V_i`, thus getting the number density :math:`\rho_{n,i}` in these
sub volumes. Inside the pore this is done by creating a radial slicing,
similar to the radial distribution function. These sub volumes are calculated by

.. math::

    V_i^\text{radial}=\pi z_\text{pore}(r_i^2-r_{i-1}^2).

with pore length :math:`z_\text{pore}` and radius :math:`r_i` of sub volume
:math:`i`. This yields

.. math::

    \rho_{n,i}^\text{radial}=\frac{N_i}{V_i^\text{radial}}=\frac{N_i}{\pi z_\text{pore}}\frac{1}{r_i^2-r_{i-1}^2}.

Outside the pore, the sub volumes are given by

.. math::

    V_j^\text{out}=(x_\text{pore}\cdot y_\text{pore}-\pi r^2)z_j

with pore width :math:`x_\text{pore}`, height :math:`y_\text{pore}`, pore radius
:math:`r` and slice width :math:`z_j`. Thus

.. math::

    \rho_{n,j}^\text{out}=\frac{N_j}{V_j^\text{out}}=\frac{N_j}{x_\text{pore}\cdot y_\text{pore}-\pi r^2}\frac{1}{z_j}.

Note that the outside refers to the reservoirs of the pore simulation.
Therefore the slices add up to the reservoir length :math:`z_{res}`.
Since there is a reservoir on each side, they are brought together
by translating the atom coordinates to one of the reservoirs. Since the
outside density refers to the density of the outside surface, it does
not contain the cylindrical extension of the pore inside the reservoirs.

Finally the mass density is calculated by

.. math::

    \rho=\frac M{N_A}\rho_n

with Avogadro constant :math:`N_A`. The units are then transformed to
:math:`\frac{\text{kg}}{\text m^3}` by

.. math::

    [\rho]=\frac{[M]\frac{\text{g}}{\text{mol}}}{[N_A]10^{23}\frac{\#}{\text{mol}}}[\rho_n]\frac{\#}{\text{nm}^3}
           =\frac{[M]}{[N_A]}[\rho_n]\cdot10\frac{\text{kg}}{\text m^3}

where the square brackets mean, that only the variables value is taken.
Since finding full molecules in a sub volume is difficult, the atoms
of the specified molecule are counted in the sub volumes and the result
is then divided by the number of atoms the molecule consists of.

.. note::

  Necessary parameters like reservoir length and pore diameter can be imported
  from the backed-up pore object.


.. raw:: html

        </div>
      </div>
    </div>
  </div>
