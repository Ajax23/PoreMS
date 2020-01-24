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

The pore generator is a class, inheriting all functionalities from the mentioned molecule builder. Accelerating optimization logic has been separated into custom classes in order to reduce clutter and simplify readability.

Finally, a storage class was implemented which translates the generated data construct into the most popular structure file-formats to be used in molecular simulations.


Molecule class
--------------

Initially the idea behind this class was an easy construction of geometries, which during development got more tuned towards molecules, containing specific information like mass and charge. Atoms can be added using specified positions or based on the location of other atoms, hereby using spherical coordinates. The underlying data structure is a matrix that contains positional information of atoms and their type.

Be :math:`\boldsymbol{D}\in\mathbb{R}^{m\times n}\mathbb{S}^{n}` a molecule with spatial dimensions :math:`i=1,\dots,m` and atom number :math:`j=1,\dots,n`. The data matrix then has following structure

.. math::

  \boldsymbol{D}=\begin{bmatrix}
  \begin{pmatrix}d_{11}\\d_{12}\\\vdots\\d_{1n}\end{pmatrix}&
  \begin{pmatrix}d_{21}\\d_{22}\\\vdots\\d_{2n}\end{pmatrix}&
  \dots&
  \begin{pmatrix}d_{m1}\\d_{m2}\\\vdots\\d_{mn}\end{pmatrix}&
  \begin{pmatrix}t_1\\t_2\\\vdots\\t_n\end{pmatrix}
  \end{bmatrix}

with position :math:`d_{ij}` of atom :math:`j` and chemical atom type :math:`t_j`. Based on these positional information, geometrical functions were implemented assuring that any desired arrangement of atoms and molecules can be realized.


Pore Class
----------

In order to create an ideal pore structure with the desired dimensions and diameter, multiple steps are involved. First, using the constructional capabilities of the molecule class, an ideal minimal :math:`\beta`-cristobalite crystal structure \ref{fig:11} was generated using the bond angle and length

.. math::

  \alpha_\text{O-Si-O}=\frac{2\cdot180}{\pi}\tan^{-1}\sqrt2=109.47^\circ,\ \ \ l_\text{Si-O}=0.155\ \text{nm}.

This is necessary to assure the requested system size is achieved by multiplying this minimal structure in all dimensions. The resulting :math:`\beta`-cristobalite block has bond lengths and angles between all silicon and oxygen atoms that are identical, and is thus feasible to be called an ideal crystal. Out of this crystal, a pore is carved out by removing all atoms within a hypothetical cylinder with the desired diameter placed inside this block.

.. figure::  /pics/struct.svg
  :align: center
  :width: 70%
  :name: fig1
  :alt: Minimal structure viewed from different perspectives.

**Figure 1:** Minimal structure viewed (a) from the :math:`xy`-plane, (b) from the :math:`yz`-plane (c) and from the :math:`xz`-plane. Silicon atoms are coloured yellow and oxygen atoms orange.

For the purpose of ensuring chemical propriety, the carved out surface needs to be processed based on a set of rules as proposed by Coasne et al. \cite{Coasne:2008}. First, all unsaturated silicon atoms are to be removed. Next, silicon atoms with three unsaturated oxygen bonds are eliminated. Finally, now unbound oxygen atoms must be deleted. The resulting surface has fully saturated silicon atoms with a maximum of two unsaturated oxygen bonds, the latter will be used as binding sites to connect molecules that are to be placed on the surface.

However, each processing step involves a high computational effort. In every step, the number of bonds for all atoms must be determined by comparing the distances of each atom to all other atoms. Typically, this effort scales quadratically with the number of atoms

.. math::

  \mathcal O(n^2).

An effort, that is adverse for highly scalable systems. The optimization itself will be discussed in the subsequent sections. A convenient by-product of the implemented optimization is the ability to carve out the cristobalite crystal in every orientational axis with an arbitrary pore shape not being limited to a cylindrical shape.

Once the surface is prepared and binding sites are determined, it is possible to attach molecules on the surface using implemented functions that automatically rotate the groups towards the central axis of the pore on the inside and orient them perpendicular to the surface on the outside. The molecules must be defined using said molecule builder.

At the end the system is finalized by saturating empty bonds with hydrogen, creating silanol groups, calculating the excess charge and positionally centring the pore system.

Store class
-----------

Molecules and pores generated using the described classes, can be converted to a readable structure in the most generic formats using the storing functionalities. These functions were separated into an own class to allow an accessible way to extend the capabilities to other needed formats.

Verlet class
------------

As mentioned before the effort for determining the bonds of each atom is unpractical, since it scales quadratically. An excellent solution is provided by the algorithms used in molecular simulations when considering short ranged interactions. These have a converging error when introducing a cut-off radius to the potential function. Atoms farther than the cut-off radius are neglected from the calculation. Similarly, bond lengths are constant. Therefore, it is only necessary to search for bonds in the close vicinity of an atom.

The implemented algorithm splits the whole system into smaller cubes of the same size which contain intersecting atoms. A scan for the bonds of an atom is then performed solely in the cube that contains the atom and the immediate neighbouring cubes, a total of 27 cubes. The computational effort for each atom is thus a constant

.. math::

  \mathcal{O}(27\cdot b^2)

since due to the ideal nature of the crystal, the number of atoms :math:`b` in a cube is constant. Therefore, the computational effort for an entire search scales linear with the number of cubes. For example, doubling the cristobalite block size only increases the effort eightfold.

Furthermore, the search is easily parallelizable, since no communication is needed between the subprocesses that each cover a set of cubes. The effort therefore has an ideal speedup.

Bonding class
-------------

Although the search can be parallelized, still multiple iterations are needed to cover the surface preparations. Additionally, due to machine inaccuracy there is the risk of bonds not being detected as such, leading to artefacts, and it is also not possible to ensure that all bonds were found when deleting atoms because all systems are shaped differently. Therefore, another optimization, or rather supporting algorithm, was implemented to bypass these issues.

The idea was reducing the number of iterations to a single search by creating a connectivity matrix of all oxygen and silicon atoms. Therefore two matrices :math:`\boldsymbol{O}` for oxygen and :math:`\boldsymbol{S}` for silicon were defined

.. math::

  \begin{array}{cc}
    \boldsymbol{O}=
      \begin{bmatrix}
          o_0&\begin{pmatrix}p_{s,0,0}&p_{s,0,1}\end{pmatrix}\\
          o_1&\begin{pmatrix}p_{s,1,0}&p_{s,1,1}\end{pmatrix}\\
          \vdots&\vdots\\
          o_k&\begin{pmatrix}p_{s,k,0}&p_{s,k,1}\end{pmatrix}
      \end{bmatrix}
    ,&
    \boldsymbol{S}=
      \begin{bmatrix}
          s_0&\begin{pmatrix}p_{o,0,0}&p_{o,0,1}&p_{o,0,2}&p_{o,0,3}\end{pmatrix}\\
          s_1&\begin{pmatrix}p_{o,1,0}&p_{o,1,1}&p_{o,1,2}&p_{o,1,3}\end{pmatrix}\\
          \vdots&\vdots\\
          s_l&\begin{pmatrix}p_{o,l,0}&p_{o,l,1}&p_{o,l,2}&p_{o,l,3}\end{pmatrix}
      \end{bmatrix}
  \end{array}


with atom ids of oxygen :math:`o` and silicon :math:`s` in the data matrix of the pore, list id pointer :math:`p_s` of the silicon entry in matrix :math:`\boldsymbol{S}` and pointer :math:`p_o` of the oxygen entry in matrix :math:`\boldsymbol{O}`. Thus each entry of the matrix presents the atom and all its binding partners. These matrices are filled after creating the cristobalite block. This way it is possible to check whether all bonds were found, since all entries need to have the same size when considering periodic boundary conditions.

Using this implementation, it is no longer required to physically delete atoms when carving out the pore, it is enough to remove binding partners from the matrices. Thus, the surface preparation conditions only need to consider the number of bonds remaining in each entry and thereby determine whether an atom needs to be removed or not, resulting into an effort scaling linear with the number of atoms

.. math::

  \mathcal{O}(n).



.. raw:: html

        </div>
      </div>
    </div>
  </div>
