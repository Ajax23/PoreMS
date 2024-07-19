################################################################################
# Molecule Class                                                               #
#                                                                              #
"""All necessary function for creating and editing molecules."""
################################################################################


import math
import pandas as pd

import porems.utils as utils
import porems.database as db
import porems.geometry as geometry

from porems.atom import Atom


class Molecule:
    """This class defines a molecule object, which is basically a list of Atom
    objects. Each atom object has a specific cartesian position and atom type,
    creating the construct of a molecule.

    Functions have been provided for editing, moving and transforming the
    atom objects, either as a collective or specific part of the molecule.

    Parameters
    ----------
    name : string, optional
        Molecule name
    short : string, optional
        Molecule short name
    inp : None, string, list, optional
        None for empty Molecule, string to read a molecule from a
        specified file link or a list of either molecules to concatenate these
        into one object, or a list of atom objects

    Examples
    --------
    Following example generates a benzene molecule without hydrogen atoms

    .. code-block:: python

        import porems as pms

        mol = pms.Molecule("benzene", "BEN")
        mol.add("C", [0,0,0])
        mol.add("C", 0, r=0.1375, theta= 60)
        mol.add("C", 1, r=0.1375, theta=120)
        mol.add("C", 2, r=0.1375, theta=180)
        mol.add("C", 3, r=0.1375, theta=240)
        mol.add("C", 4, r=0.1375, theta=300)

    """
    def __init__(self, name="molecule", short="MOL", inp=None):
        # Initialize
        self._dim = 3

        self._name = name
        self._short = short

        self._box = []
        self._charge = 0
        self._masses = []
        self._mass = 0

        # Check data input
        if inp is None:
            self._atom_list = []
        else:
            # Read from file
            if isinstance(inp, str):
                self._atom_list = self._read(inp, inp.split(".")[-1].upper())
            # Concat multiple molecules
            elif isinstance(inp, list):
                # Atom list is provided
                if(isinstance(inp[0], Atom)):
                    self._atom_list = inp
                # List of molecules is provided
                else:
                    self._atom_list = self._concat(inp)


    ##################
    # Representation #
    ##################
    def __repr__(self):
        """Create a pandas table of the molecule data.

        Returns
        -------
        repr : DataFrame
            Pandas data frame of the molecule object
        """
        # Set colums names
        columns = ["Residue", "Name", "Type", "x", "y", "z"]

        # Get data
        data =[]
        for atom in self._atom_list:
            data.append([atom.get_residue(), atom.get_name(), atom.get_atom_type(),
                         atom.get_pos()[0], atom.get_pos()[1], atom.get_pos()[2]])

        # Create data frame
        return pd.DataFrame(data, columns=columns).to_string()


    ##############
    # Management #
    ##############
    def _read(self, file_path, file_type):
        """Read a molecule from a file. Currently only **GRO**, **PDB** and
        **MOL2** files are supported.

        Parameters
        ----------
        file_path : string
            Link to requested file
        file_type : string
            File extension name

        Returns
        -------
        atom_list : list
            Atom list
        """
        # Process input
        if not file_type in ["GRO", "PDB", "MOL2"]:
            print("Unsupported filetype.")
            return

        # Read molecule
        atom_list = []
        with open(file_path, "r") as file_in:
            for line_idx, line in enumerate(file_in):
                line_val = line.split()
                is_add = False

                # Gro file
                if file_type == "GRO":
                    if line_idx > 0 and len(line_val) > 3:
                        residue = int(line[0:5])-1
                        pos = [float(line_val[i]) for i in range(3, 5+1)]
                        name = line_val[1]
                        atom_type = ''.join([i for i in line_val[1] if not i.isdigit()])
                        is_add = True

                # Pdb file
                elif file_type == "PDB":
                    if line_val[0] in ["ATOM", "HETATM"]:
                        residue = int(line[22:26])-1
                        pos = [float(line_val[i])/10 for i in range(6, 8+1)]
                        name = line_val[11]
                        atom_type = line_val[11]
                        is_add = True

                # Mol2 file
                elif file_type == "MOL2":
                    if len(line_val) > 8:
                        residue = 0
                        pos = [float(line_val[i])/10 for i in range(2, 4+1)]
                        name = line_val[1]
                        atom_type = ''.join([i for i in line_val[1] if not i.isdigit()])
                        is_add = True

                if is_add:
                    atom_list.append(Atom(pos, atom_type, name, residue))

        # Transform to column
        return atom_list

    def _concat(self, mol_list):
        """Concatenate a molecule list into one molecule object.

        Parameters
        ----------
        mol_list : list
            List of molecule objects to be concatenated

        Returns
        -------
        atom_list : list
            Atom list
        """
        return sum([mol.get_atom_list() for mol in mol_list], [])

    def _temp(self, atoms):
        """Create a temporary molecule of specified atom ids.

        Parameters
        ----------
        atoms : list
            List of atoms to be included

        Returns
        -------
        mol : Molecule
            Molecule object
        """
        return Molecule(inp=[self._atom_list[x] for x in atoms])

    def append(self, mol):
        """Append a given molecule to the current object.

        Parameters
        ----------
        mol : Molecule
            Molecule object
        """
        self._atom_list += mol.get_atom_list()

    def column_pos(self):
        """Create column list of atom positions

        Returns
        -------
        column : list
            Columns of all atom positions in all dimensions
        """
        return utils.column([atom.get_pos() for atom in self._atom_list])


    ############
    # Geometry #
    ############
    def _vector(self, pos_a, pos_b):
        """Calculate the vector between to two positions as defined in
        :class:`porems.geometry.vector` with the addition to define the inputs
        as atom indices.

        Parameters
        ----------
        pos_a : integer, list
            First position :math:`\\boldsymbol{a}`
        pos_b : integer, list
            Second position :math:`\\boldsymbol{b}`

        Returns
        -------
        vector : list
            Bond vector
        """
        # Process input
        if isinstance(pos_a, int) and isinstance(pos_b, int):
            pos_a = self.pos(pos_a)
            pos_b = self.pos(pos_b)
        elif not (isinstance(pos_a, list) and isinstance(pos_b, list)):
            print("Vector: Wrong input...")
            return

        # Check dimensions
        if not len(pos_a) == self._dim:
            print("Vector: Wrong dimensions...")
            return

        # Calculate vector
        return geometry.vector(pos_a, pos_b)

    def _box_size(self):
        """Calculate the box size of the current molecule. This is done by
        determining the maximal coordinate value of all atoms in all dimensions

        .. math::

            \\boldsymbol{b}=\\begin{pmatrix}\\max(\\boldsymbol{d}_1)&\\max(\\boldsymbol{d}_1)&\\dots&\\max(\\boldsymbol{d}_n)\\end{pmatrix}^T

        where :math:`\\boldsymbol{d}_i` is the dimension-vector of the data
        matrix.

        Returns
        -------
        box : list
            Box length of the current molecule
        """
        data = self.column_pos()
        return [max(data[i]) if max(data[i]) > 0 else 0.001 for i in range(self._dim)]


    ##############
    # Properties #
    ##############
    def pos(self, atom):
        """Get the position of an atom.

        Parameters
        ----------
        atom : integer
            Atom id

        Returns
        -------
        pos : list
            Position vector of the specified atom
        """
        return self._atom_list[atom].get_pos()

    def bond(self, inp_a, inp_b):
        """Return the bond vector of a specified bond. The two inputs can either
        be atom indices or to vectoral positions.

        Parameters
        ----------
        inp_a : integer, list
            Either an atom id or a position vector
        inp_b : integer, list
            Either an atom id or a position vector

        Returns
        -------
        bond : list
            Bond vector

        Examples
        --------
        .. code-block:: python

            mol.bond(0, 1)
            mol.bond(*[0, 1])
            mol.bond([1, 0, 0], [0, 0, 0])
        """
        return self._vector(inp_a, inp_b)

    def centroid(self):
        """Calculate the geometrical center of mass

        .. math::

            \\text{centroid}=\\begin{pmatrix}c_1&c_2&\\dots&c_n\\end{pmatrix}^T

        with

        .. math::

            c_i=\\frac{1}{m}\\sum_j^m d_{ij}.

        Hereby :math:`i\\dots n` stands for the dimension and :math:`j\\dots m`
        for the molecule.

        Returns
        -------
        Centroid : list
            Geometrical center of mass
        """
        # Calculate the centroid
        data = self.column_pos()
        return [sum(data[i])/len(data[i]) for i in range(self._dim)]

    def com(self):
        """Calculate the center of mass

        .. math::

            \\text{com}=\\begin{pmatrix}c_1&c_2&\\dots&c_n\\end{pmatrix}^T

        with

        .. math::

            c_i=\\frac{1}{\\sum_j^mM_j}\\sum_j^m d_{ij}\\cdot M_j

        and mass :math:`M`. Hereby :math:`i\\dots n` stands for the dimension
        and :math:`j\\dots m` for the molecule.

        Returns
        -------
        COM : list
            Center of mass
        """
        # Calculate the center of mass
        data = self.column_pos()
        masses = self.get_masses()
        return [sum([data[i][j]*masses[j] for j in range(self.get_num())])/sum(masses) for i in range(self._dim)]


    #################
    # Basic Editing #
    #################
    def translate(self, vec):
        """Translate the atoms data matrix :math:`\\boldsymbol{D}` along a
        vector :math:`\\boldsymbol{a}\\in\\mathbb{R}^n`.

        .. math::

            \\boldsymbol{D}_\\text{trans}=
            \\boldsymbol{D}+\\boldsymbol{a}=
            \\begin{pmatrix}
            \\boldsymbol{d}_1+a_1&\\boldsymbol{d}_2+a_2&\\dots&\\boldsymbol{d}_n+a_n&\\boldsymbol{d}_t
            \\end{pmatrix}

        Parameters
        ----------
        vec : list
            Vector a
        """
        for atom in self._atom_list:
            atom.set_pos([atom.get_pos()[i]+vec[i] for i in range(self._dim)])

    def rotate(self, axis, angle, is_deg=True):
        """Rotate data matrix :math:`\\boldsymbol{D}` around an axis
        :math:`\\boldsymbol{a}\\in\\mathbb{R}^3` with angle
        :math:`\\alpha\\in\\mathbb{R}` using rotation function
        :func:`porems.geometry.rotate`.

        Parameters
        ----------
        axis : integer, string, list
            Axis
        angle : float
            Angle
        is_deg : bool, optional
            True if the input is in degree
        """
        for atom in self._atom_list:
            atom.set_pos(geometry.rotate(atom.get_pos(), axis, angle, is_deg))

    def move(self, atom, pos):
        """Move whole the molecule to a new position, where the dragging point
        is a given atom that is moved to a specified position.

        Parameters
        ----------
        atom : integer
            Main atom id whose position will be changed
        pos : list
            New position vector
        """
        self.translate(self._vector(self.pos(atom), pos))

    def zero(self, pos=[0, 0, 0]):
        """Move whole the molecule, so that the minimal coordinate
        between all atoms is zero in all dimensions, or rather the values of the
        position variable ``pos``. This function is basically setting
        the zero point of the coordinate system to ``pos``.

        Parameters
        ----------
        pos : list, optional
            Vector of the zero point of the coordinate system

        Returns
        -------
        vec : list
            Vector used for the translation
        """
        # Calculate translation vector
        data = self.column_pos()
        vec = [pos[i]-min(data[i]) for i in range(self._dim)]

        # Reset box size
        self._box = []

        # Translate molecule
        self.translate(vec)

        return vec

    def put(self, atom, pos):
        """Change the position of an atom to a given position vector
        :math:`\\boldsymbol{a}\\in\\mathbb{R}^n`.

        Parameters
        ----------
        atom : integer
            Atom id whose position will be changed
        pos : list
            New position vector
        """
        self._atom_list[atom].set_pos(pos)


    ####################
    # Advanced Editing #
    ####################
    def part_move(self, bond, atoms, length, vec=[]):
        """Change the length of a specified bond. Variable ``atoms`` specifies
        which atoms or rather which part of the molecule needs to be moved for
        this specific bond. The given length is going to be the new bond length,
        **not** by how much the bond length is changed.

        The move vector is determined automatically by the given length and atom
        bond. This vector can also be given manually with no regards to length,
        by setting the variable ``vec``.

        Parameters
        ----------
        bond : list
            List of two atom ids of the bond to be adjusted
        atoms : integer, list
            List of atoms that need to be moved by changing the bond length
            (can also be one id)
        length : float
            New bond length
        vec : list, optional
            Set this vector to manually define the translation vector

        Examples
        --------
        .. code-block:: python

            mol.part_move([0, 1], [1, 2, 3], 0.5)
        """
        # Create temporary molecule
        if isinstance(atoms, int):
            atoms = [atoms]
        temp = self._temp(atoms)

        # Set length
        length = abs(length-geometry.length(self.bond(*bond)))

        # Set vector
        if not vec:
            vec = self._vector(bond[0], bond[1])
        vec = [v*length for v in geometry.unit(vec)]

        # Move molecule
        temp.translate(vec)

    def part_rotate(self, bond, atoms, angle, zero):
        """Rotate a set of specified atoms around a given bond as the rotation
        axis. First however the system needs to be set to zero. Therefore the
        atom id to define the new coordinate system must be given for the set
        of specified atoms. Normally this is the atoms that is on the end of the
        given bond axis.

        Parameters
        ----------
        bond : list
            List of two atom ids of the bond to be set as an axis
        atoms : integer, list
            List of atoms to be rotated (can also be one id)
        angle : float
            Rotation angle
        zero : integer
            Atom id to define zero point of the new coordinate system

        Examples
        --------
        .. code-block:: python

            mol.part_rotate([0, 1], [1, 2, 3], 90, 0)
        """
        # Create temporary molecule
        self.move(zero, [0, 0, 0])
        if isinstance(atoms, int):
            atoms = [atoms]
        temp = self._temp(atoms)

        # Rotate molecule
        temp.rotate([self.pos(bond[0]), self.pos(bond[1])], angle)

    def part_angle(self, bond_a, bond_b, atoms, angle, zero):
        """Change the bond angle of two bond vectors. Variable ``atoms``
        specifies which atoms or rather which part of the molecule needs to be
        rotated in order to change the specified bond angle. First however the
        system needs to be set to zero. Therefore, the atom id to define the new
        coordinate system must be given for the set of specified atoms.
        Normally this is the atom that touches the angle.

        The rotation axis is determined by creating the cross product
        of the two bond vectors. Thus, getting the normal vector of a surface
        that contains both bond vectors.

        Parameters
        ----------
        bond_a : list
            First bond vector given as a list of two atom ids
        bond_b : list
            Second bond vector given as a list of two atom ids
        atoms : integer, list
            List of atoms to be rotated (can also be one id)
        angle : float
            Rotation angle
        zero : integer
            Atom id to define zero point of the new coordinate system

        Examples
        --------
        .. code-block:: python

            mol.part_angle([0, 1], [1, 2], [1, 2, 3], 90, 1)
        """
        # Create temporary molecule
        self.move(zero, [0, 0, 0])
        if isinstance(atoms, int):
            atoms = [atoms]
        temp = self._temp(atoms)

        # Rotate molecule around normal vector
        if len(bond_a) == len(bond_b):
            if len(bond_a) == 2:
                vec = geometry.cross_product(self._vector(*bond_a), self._vector(*bond_b))
            elif len(bond_a) == self._dim:
                vec = geometry.cross_product(bond_a, bond_b)
            else:
                print("Part_Angle : Wrong bond input...")
                return
        else:
            print("Part_Angle : Wrong bond dimensions...")
            return

        # Rotate molecule
        temp.rotate(vec, angle)


    #########
    # Atoms #
    #########
    def add(self, atom_type, pos, bond=[], r=0, theta=0, phi=0, is_deg=True, name="", residue=0):
        """Add a new atom in polar coordinates. The ``pos`` input is either
        an atom id, that determines is the bond-start,
        or a vector for a specific position.

        If the polar coordinates are dependent on a bond vector as an axis,
        the ``bond`` variable must be set. The coordinate system is then
        transformed to the bond axis. If this variable is set to an empty list,
        then the given coordinates are assumed to be dependent on the z-axis.

        Parameters
        ----------
        atom_type : string
            Atom type
        pos : integer, list
            Position of the atom
        bond : list, optional
            Bond axis
        r : float, optional
            Bond length
        theta : float, optional
            Azimuthal angle
        phi : float, optional
            Polar angle
        is_deg : bool, optional
            True if the input of the angles in degree
        name : string, optional
            Optionally set a unique atom name
        residue : integer, optional
            Residue number

        Examples
        --------
        .. code-block:: python

            mol.add("C", [0, 0, 0])
            mol.add("C", 0, r=0.153, theta=-135)
            mol.add("C", 1, [0, 1], r=0.153, theta= 135)
        """
        # Process input
        pos = self.pos(pos) if isinstance(pos, int) else pos
        vec = self._vector(*bond) if bond else geometry.main_axis("z")

        # Add coordinate transformation when given a bond
        phi += geometry.angle_polar(vec, is_deg)
        theta += geometry.angle_azi(vec, is_deg)

        # Process angles
        phi *= math.pi/180 if is_deg else 1
        theta *= math.pi/180 if is_deg else 1

        # Transform spherical to cartesian coordinates
        x = r*math.sin(theta)*math.cos(phi)
        y = r*math.sin(theta)*math.sin(phi)
        z = r*math.cos(theta)
        coord = [x, y, z]

        # Create new atom
        self._atom_list.append(Atom([pos[i]+coord[i] for i in range(self._dim)], atom_type, name, residue))

    # Delete an atom
    def delete(self, atoms):
        """Delete specified atom from the molecule. The input can also be a list
        of atom ids.

        Parameters
        ----------
        atoms : integer, list
            Atom id or list to be deleted
        """
        # Process input
        atoms = [atoms] if isinstance(atoms, int) else atoms

        # Remove atoms
        for atom in sorted(atoms, reverse=True):
            self._atom_list.pop(atom)

    def overlap(self, error=0.005):
        """Search for overlapping atoms.

        Parameters
        ----------
        error : float, optional
            Error tolerance

        Returns
        -------
        duplicate : dictionary
            Dictionary of duplicate lists
        """
        # Initialize
        atom_list = {x: False for x in range(self.get_num())}
        duplicates = {}

        # Run through complete atoms list
        for atom_a in atom_list:
            # Ignore duplicate items
            if not atom_list[atom_a]:
                # Run through atom list after first loop
                for atom_b in [x for x in atom_list if x>atom_a]:
                    # Check if overlapping
                    if sum([error > abs(x) for x in geometry.vector(self.pos(atom_a), self.pos(atom_b))]) == 3:
                        if not atom_a in duplicates:
                            duplicates[atom_a] = []
                        duplicates[atom_a].append(atom_b)
                        # Set to false
                        atom_list[atom_a] = True
                        atom_list[atom_b] = True

        # Return duplicates
        return duplicates

    def switch_atom_order(self, atom_a, atom_b):
        """Switch atom order of two atoms.

        Parameters
        ----------
        atom_a : integer
            Atom id of the first atom
        atom_b : integer
            Atom id of the second atom
        """
        self._atom_list[atom_a], self._atom_list[atom_b] = self._atom_list[atom_b], self._atom_list[atom_a]

    def set_atom_type(self, atom, atom_type):
        """Change the atom type of a specified atom.

        Parameters
        ----------
        atom : integer
            Atom id
        atom_type : string
            New atom type
        """
        self._atom_list[atom].set_atom_type(atom_type)

    def set_atom_name(self, atom, name):
        """Change the atom name of a specified atom.

        Parameters
        ----------
        atom : integer
            Atom id
        name : string
            New atom name
        """
        self._atom_list[atom].set_name(name)

    def set_atom_residue(self, atom, residue):
        """Change the residue index of a specified atom.

        Parameters
        ----------
        atom : integer
            Atom id
        residue : integer
            New residue index
        """
        self._atom_list[atom].set_residue(residue)

    def get_atom_type(self, atom):
        """Return the atom type of the given atom id.

        Parameters
        ----------
        atom : integer
            Atom id

        Returns
        -------
        atom_type : string
            Atom type
        """
        return self._atom_list[atom].get_atom_type()

    def get_atom_list(self):
        """Return the atoms list.

        Returns
        -------
        atom_list : list
            Atom list
        """
        return self._atom_list


    ##################
    # Setter Methods #
    ##################
    def set_name(self, name):
        """Set the molecule name.

        Parameters
        ----------
        name : string
            Molecule name
        """
        self._name = name

    def set_short(self, short):
        """Set the molecule short name.

        Parameters
        ----------
        short : string
            Molecule short name
        """
        self._short = short

    def set_box(self, box):
        """Set the box size.

        Parameters
        ----------
        box : list
            Box size in all dimensions
        """
        self._box = box

    def set_charge(self, charge):
        """Set the total charge of the molecule.

        Parameters
        ----------
        charge : float
            Total molecule charge
        """
        self._charge = charge

    def set_masses(self, masses=[]):
        """Set the molar masses of the atoms.

        Parameters
        ----------
        masses : list, optional
            List of molar masses in :math:`\\frac g{mol}`
        """
        self._masses = masses if masses else [db.get_mass(atom.get_atom_type()) for atom in self._atom_list]

    def set_mass(self, mass=0):
        """Set the molar mass of the molecule.

        Parameters
        ----------
        mass : float, optional
            Molar mass in :math:`\\frac g{mol}`
        """
        self._mass = mass if mass else sum(self.get_masses())


    ##################
    # Getter Methods #
    ##################
    def get_name(self):
        """Return the molecule name.

        Returns
        -------
        name : string
            Molecule name
        """
        return self._name

    def get_short(self):
        """Return the molecule short name.

        Returns
        -------
        short : string
            Molecule short name
        """
        return self._short

    def get_box(self):
        """Return the box size of the molecule.

        Returns
        -------
        box : list
            Box size in all dimensions
        """
        return self._box if self._box else self._box_size()

    def get_num(self):
        """Return the number of atoms.

        Returns
        -------
        num : integer
            Number of atoms
        """
        return len(self._atom_list)

    def get_charge(self):
        """Return the total charge of the molecule.

        Returns
        -------
        charge : float
            Total charge
        """
        return self._charge

    def get_masses(self):
        """Return a list of masses of the atoms.

        Returns
        -------
        masses : list
            Masses in :math:`\\frac g{mol}`
        """
        if not self._masses:
            self.set_masses()
        return self._masses

    def get_mass(self):
        """Return the molar mass of the molecule.

        Returns
        -------
        mass : float
            Molar mass in :math:`\\frac g{mol}`
        """
        if not self._mass:
            self.set_mass()
        return self._mass
