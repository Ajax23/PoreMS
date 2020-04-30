:orphan:

.. raw:: html

  <div class="container-fluid">
    <div class="row">
      <div class="col-md-10">
        <div style="text-align: justify; text-justify: inter-word;">

Operating principle
===================

A pore generation package for molecular simulations (PoreMS) has been programmed in an object-oriented manner using python 3. The main idea was maximizing individual customizational possibilities without obscuring usability.
For this purpose, the core piece of this tool is a molecule construction set containing functions to place, translate, and rotate atoms in a three-dimensional space. Hence custom molecules and systems can be created in the required specification, for example moving a molecule into the binding pocket of an enzyme as a starting configuration of a simulation.

The construction of a pore takes place in multiple steps from creating a grid pattern, carving out geometrical shapes and functionalize the surface. Accelerating optimization logic has been separated into custom classes in order to reduce clutter and simplify readability.

Finally, a storage class was implemented which translates the generated data construct into the most popular structure file-formats to be used in molecular simulations.


Molecule Generation
-------------------

Initially the idea behind the molecule class was an easy construction of geometries, which during development got more tuned towards molecules, containing specific information like mass and charge. Atoms can be added using specified positions or based on the location of other atoms, hereby using spherical coordinates. The underlying data structure is a list of atom objects, with each of these objects containing the atom position, atom type, and atom name. Based on these positional information, geometrical functions were implemented assuring that any desired arrangement of atoms and molecules can be realized.


Pore Generation
---------------

In order to create an ideal pore structure with the desired dimensions and diameter, multiple steps are involved.


Pattern
~~~~~~~

First, using the constructional capabilities of the molecule class, a pattern class has been devised, which creates a minimal :math:`\beta`-cristobalite crystal structure \ref{fig:11} using the bond angle and length

.. math::

  \alpha_\text{O-Si-O}=\frac{2\cdot180}{\pi}\tan^{-1}\sqrt2=109.47^\circ,\ \ \ l_\text{Si-O}=0.155\ \text{nm}.

This is necessary to assure the requested system size is achieved as best as possible by multiplying this minimal structure in all dimensions. The resulting :math:`\beta`-cristobalite block has bond lengths and angles between all silicon and oxygen atoms that are identical, and is thus feasible to be called an ideal crystal. Using the provided pattern class as a basis, other structural patterns can be created if needed.

.. figure::  /pics/struct.svg
  :align: center
  :width: 70%
  :name: fig1
  :alt: Minimal structure viewed from different perspectives.

**Figure 1:** Minimal structure of :math:`\beta`-cristobalite viewed (a) from the :math:`xy`-plane, (b) from the :math:`yz`-plane (c) and from the :math:`xz`-plane. Silicon atoms are colored yellow and oxygen atoms orange.


Shape
~~~~~

Out of this molecule block, a volume is carved out by removing all atoms within a hypothetical shape placed inside the system. Sub classes of the pattern class have been created defining the basic geometrical forms like cylindrical and spherical shapes. Other volumetric shapes can be created by combining provided classes or by defining own shapes using the provided ones as an example.


Preparation
~~~~~~~~~~~

For the purpose of ensuring chemical propriety, the carved out surface needs to be processed based on a set of rules as proposed by Coasne et al. \cite{Coasne:2008}. First, all unsaturated silicon atoms are to be removed. Next, silicon atoms with three unsaturated oxygen bonds are eliminated. Finally, now unbound oxygen atoms must be deleted. The resulting surface has fully saturated silicon atoms with a maximum of two unsaturated oxygen bonds, the latter will be used as binding sites to connect molecules that are to be placed on the surface.

However, each processing step involves a high computational effort. In every step, the number of bonds for all atoms must be determined by comparing the distances of each atom to all other atoms. Typically, this effort scales quadratically with the number of atoms

.. math::

  \mathcal O(n^2).

An effort, that is adverse for highly scalable systems. Once the surface is prepared and binding sites are determined, it is possible to attach molecules on the surface using implemented functions that automatically rotate the groups perpendicular to surface using normal vectors implemented in the shape classes.


Dice
~~~~

As mentioned, the effort for determining the bonds of each atom is unpractical, since it scales quadratically. An excellent solution is provided by the algorithms used in molecular simulations when considering short ranged interactions. These have a converging error when introducing a cut-off radius to the potential function. Atoms farther than the cut-off radius are neglected from the calculation. Similarly, bond lengths are constant. Therefore, it is only necessary to search for bonds in the close vicinity of an atom.

The implemented algorithm splits the whole system into smaller cubes of the same size which contain intersecting atoms. A scan for the bonds of an atom is then performed solely in the cube that contains the atom and the immediate neighboring cubes, a total of 27 cubes. The computational effort for each atom is thus a constant

.. math::

  \mathcal{O}(27\cdot b^2)

since due to the geometrical ideal nature of the crystal, the number of atoms :math:`b` in a cube is constant. Therefore, the computational effort for an entire search scales linear with the number of cubes. For example, doubling the cristobalite block size only increases the effort eightfold.

Furthermore, the search is easily parallelizable, since no communication is needed between the subprocesses that each cover a set of cubes. The effort therefore has an ideal speedup.


Matrix
~~~~~~

Although the search can be parallelized, still multiple iterations are needed to cover the surface preparations. Additionally, due to machine inaccuracy there is the risk of bonds not being detected as such, leading to artefacts. Also, it is not possible to ensure that all bonds were found when deleting atoms, because all systems are shaped differently. Therefore, another optimization, or rather supporting algorithm, was implemented to bypass these issues.

The idea was reducing the number of iterations to a single search by creating a connectivity matrix of all grid atoms. The result is a dictionary that has atoms :math:`1\dots n` as keys and their corresponding value is a list of bonded atoms :math:`1\dots m`

.. math::

  \boldsymbol{C}=
  \begin{Bmatrix}
      a_1:&\begin{bmatrix}a_{1,1}&a_{1,2}&\dots&a_{1,m_1}\end{bmatrix}\\
      a_2:&\begin{bmatrix}a_{2,1}&a_{2,2}&\dots&a_{2,m_2}\end{bmatrix}\\
      \vdots&\vdots\\
      a_n:&\begin{bmatrix}a_{n,1}&a_{n,2}&\dots&a_{n,m_n}\end{bmatrix}\\
  \end{Bmatrix}

Using this implementation, it is no longer required to physically delete atoms when carving out a structure, it is enough to remove binding partners from the matrix. For example, conditions for the surface preparation only need to consider the number of bonds remaining in each entry and thereby determine whether an atom needs to be removed or not, resulting into a negligible computational effort scaling linear with the number of atoms

.. math::

  \mathcal{O}(n).

A convenient by-product of the implemented optimizations is the ability to carve out the cristobalite crystal in every orientational axis with an arbitrary pore shape.


Store
-----

Molecules and pores generated using the described classes, can be converted to a readable structure in the most generic formats using the storing functionalities. These functions were separated into an own class to allow an accessible way to extend the capabilities to other needed formats.


.. raw:: html

        </div>
      </div>
    </div>
  </div>
