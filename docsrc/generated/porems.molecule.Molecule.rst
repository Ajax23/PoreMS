.. raw:: html

    </div>
    <div class=col-md-9 content>

Molecule
========

.. currentmodule:: porems.molecule

.. autoclass:: Molecule


  .. rubric:: Representation

  .. autosummary::

    ~Molecule.__repr__
    ~Molecule.__str__


  .. rubric:: Management

  .. autosummary::

    ~Molecule._concat
    ~Molecule._read
    ~Molecule._temp

    ~Molecule.append
    ~Molecule.column_pos


  .. rubric:: Geometry

  .. autosummary::

    ~Molecule._box_size
    ~Molecule._vector


  .. rubric:: Properties

  .. autosummary::

    ~Molecule.bond
    ~Molecule.centroid
    ~Molecule.com
    ~Molecule.pos


  .. rubric:: Editing

  .. autosummary::

    ~Molecule.move
    ~Molecule.put
    ~Molecule.rotate
    ~Molecule.translate
    ~Molecule.zero
    ~Molecule.part_angle
    ~Molecule.part_move
    ~Molecule.part_rotate


  .. rubric:: Atoms

  .. autosummary::

    ~Molecule.add
    ~Molecule.delete
    ~Molecule.overlap
    ~Molecule.switch_atom_order
    ~Molecule.set_atom_name
    ~Molecule.set_atom_type
    ~Molecule.get_atom_type


  .. rubric:: Setter Methods

  .. autosummary::

    ~Molecule.set_box
    ~Molecule.set_charge
    ~Molecule.set_mass
    ~Molecule.set_masses
    ~Molecule.set_name
    ~Molecule.set_short


  .. rubric:: Getter Methods

  .. autosummary::

    ~Molecule.get_atom_list
    ~Molecule.get_box
    ~Molecule.get_charge
    ~Molecule.get_mass
    ~Molecule.get_masses
    ~Molecule.get_max
    ~Molecule.get_name
    ~Molecule.get_num
    ~Molecule.get_short
