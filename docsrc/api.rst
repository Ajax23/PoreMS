:orphan:

API reference
=============

.. currentmodule:: porems


.. _generation_api:

Generation
----------

.. autosummary::
    :toctree: generated/

    pore.Pore
    molecule.Molecule
    molecule.Molecule.__repr__
    molecule.Molecule.__str__
    molecule.Molecule._concat
    molecule.Molecule._read
    molecule.Molecule._temp
    molecule.Molecule.append
    molecule.Molecule.column_pos
    molecule.Molecule._box_size
    molecule.Molecule._vector
    molecule.Molecule.bond
    molecule.Molecule.centroid
    molecule.Molecule.com
    molecule.Molecule.pos
    molecule.Molecule.move
    molecule.Molecule.put
    molecule.Molecule.rotate
    molecule.Molecule.translate
    molecule.Molecule.zero
    molecule.Molecule.part_angle
    molecule.Molecule.part_move
    molecule.Molecule.part_rotate
    molecule.Molecule.add
    molecule.Molecule.delete
    molecule.Molecule.overlap
    molecule.Molecule.switch_atom_order
    molecule.Molecule.set_atom_name
    molecule.Molecule.get_atom_list
    molecule.Molecule.set_atom_type
    molecule.Molecule.get_atom_type
    molecule.Molecule.set_box
    molecule.Molecule.set_charge
    molecule.Molecule.set_mass
    molecule.Molecule.set_masses
    molecule.Molecule.set_name
    molecule.Molecule.set_short
    molecule.Molecule.get_box
    molecule.Molecule.get_charge
    molecule.Molecule.get_mass
    molecule.Molecule.get_masses
    molecule.Molecule.get_name
    molecule.Molecule.get_num
    molecule.Molecule.get_short
    atom.Atom
    atom.Atom.__repr__
    atom.Atom.__str__
    atom.Atom.set_atom_type
    atom.Atom.set_name
    atom.Atom.set_pos
    atom.Atom.get_atom_type
    atom.Atom.get_name
    atom.Atom.get_pos
    store.Store
    store.Store.gro
    store.Store.pdb
    store.Store.xyz
    store.Store.obj
    store.Store.job
    store.Store.grid
    store.Store.top


.. _pattern_api:

Pattern
-------

.. autosummary::
    :toctree: generated/

    pattern.Pattern
    pattern.Pattern._block
    pattern.Pattern._orientation
    pattern.Pattern.generate
    pattern.Pattern.get_block
    pattern.Pattern.get_gap
    pattern.Pattern.get_repeat
    pattern.Pattern.get_size
    pattern.BetaCristobalit
    pattern.BetaCristobalit._hexagonal
    pattern.BetaCristobalit.pattern


.. _shape_api:

Shape
-----

.. autosummary::
    :toctree: generated/

    shape.Shape
    shape.Shape.convert
    shape.Shape.plot
    shape.Cylinder
    shape.Cylinder.Phi
    shape.Cylinder.d_Phi_phi
    shape.Cylinder.d_Phi_z
    shape.Cylinder.rim
    shape.Cylinder.surf
    shape.Cylinder.is_in
    shape.Cylinder.normal
    shape.Cylinder.surface
    shape.Cylinder.volume
    shape.Sphere
    shape.Sphere.Phi
    shape.Sphere.d_Phi_phi
    shape.Sphere.d_Phi_theta
    shape.Sphere.rim
    shape.Sphere.surf
    shape.Sphere.is_in
    shape.Sphere.normal
    shape.Sphere.surface
    shape.Sphere.volume


.. _optimization_api:

Optimization
------------

.. autosummary::
    :toctree: generated/

    dice.Dice
    dice.Dice._fill
    dice.Dice._pos_to_index
    dice.Dice._split
    dice.Dice._right
    dice.Dice._left
    dice.Dice._top
    dice.Dice._bot
    dice.Dice._front
    dice.Dice._back
    dice.Dice._step
    dice.Dice.neighbour
    dice.Dice.find_bond
    dice.Dice.find_parallel
    dice.Dice.set_pbc
    dice.Dice.get_count
    dice.Dice.get_mol
    dice.Dice.get_origin
    dice.Dice.get_pointer
    dice.Dice.get_size
    matrix.Matrix
    matrix.Matrix.bound
    matrix.Matrix.split
    matrix.Matrix.strip
    matrix.Matrix.get_matrix


.. _essential_api:

Essential Molecules
-------------------

.. autosummary::
    :toctree: generated/

    essentials.Alkane
    essentials.Alcohol
    essentials.Ketone
    essentials.TMS


.. _utils_api:

Utilities
---------

.. autosummary::
    :toctree: generated/

    utils
    utils.column
    utils.copy
    utils.load
    utils.mkdirp
    utils.replace
    utils.save
    utils.tic
    utils.toc
    geometry
    geometry.angle
    geometry.angle_azi
    geometry.angle_polar
    geometry.cross_product
    geometry.dot_product
    geometry.length
    geometry.main_axis
    geometry.rotate
    geometry.unit
    geometry.vector
    database
    database.get_mass
