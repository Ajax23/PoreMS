:orphan:

API reference
=============

.. currentmodule:: porems


.. _generation_api:

Molecule
--------

.. autosummary::
    :toctree: generated/

    atom.Atom
    atom.Atom.__repr__
    atom.Atom.__str__
    atom.Atom.set_atom_type
    atom.Atom.set_name
    atom.Atom.set_pos
    atom.Atom.get_atom_type
    atom.Atom.get_name
    atom.Atom.get_pos
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
    molecule.Molecule.part_angle
    molecule.Molecule.part_move
    molecule.Molecule.part_rotate
    molecule.Molecule.put
    molecule.Molecule.rotate
    molecule.Molecule.translate
    molecule.Molecule.zero
    molecule.Molecule.add
    molecule.Molecule.delete
    molecule.Molecule.get_atom_type
    molecule.Molecule.overlap
    molecule.Molecule.set_atom_name
    molecule.Molecule.set_atom_type
    molecule.Molecule.switch_atom_order
    molecule.Molecule.set_box
    molecule.Molecule.set_charge
    molecule.Molecule.set_mass
    molecule.Molecule.set_masses
    molecule.Molecule.set_name
    molecule.Molecule.set_short
    molecule.Molecule.get_atom_list
    molecule.Molecule.get_box
    molecule.Molecule.get_charge
    molecule.Molecule.get_mass
    molecule.Molecule.get_masses
    molecule.Molecule.get_name
    molecule.Molecule.get_num
    molecule.Molecule.get_short
    store.Store
    store.Store.obj
    store.Store.gro
    store.Store.job
    store.Store.pdb
    store.Store.xyz
    store.Store.grid
    store.Store.top


.. _pore_api:

Pore
----

.. autosummary::
    :toctree: generated/

    pore.Pore

Pattern
~~~~~~~

.. autosummary::
    :toctree: generated/

    pattern.Pattern
    pattern.Pattern.pattern
    pattern.Pattern.get_gap
    pattern.Pattern.get_repeat
    pattern.BetaCristobalit

Volume
~~~~~~

.. autosummary::
    :toctree: generated/


Optimization
~~~~~~~~~~~~

.. autosummary::
    :toctree: generated/

    bonding.Bonding
    bonding.Bonding._geminal
    bonding.Bonding._unbind
    bonding.Bonding.attach
    bonding.Bonding.delete
    bonding.Bonding.drill
    bonding.Bonding.prepare
    bonding.Bonding.remove
    bonding.Bonding.site
    bonding.Bonding.unlink
    verlet.Verlet
    verlet.Verlet._backward
    verlet.Verlet._fill
    verlet.Verlet._find_bond
    verlet.Verlet._forward
    verlet.Verlet._index
    verlet.Verlet._input
    verlet.Verlet._position
    verlet.Verlet._verlet
    verlet.Verlet.delete
    verlet.Verlet.find_parallel
    verlet.Verlet.neighbour
    verlet.Verlet.set_pbc
    verlet.Verlet.get_box
    verlet.Verlet.get_count
    verlet.Verlet.get_mol
    verlet.Verlet.get_size


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
    utils.copy_dir
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
