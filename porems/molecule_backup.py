################################################################################
# Molecule Class                                                               #
#                                                                              #
"""All necessary function for creating and editing molecules."""
################################################################################


import math

import porems.utils as utils


class Molecule:
    """This class defines a molecule object, which is basically a data matrix
    that contains atom information of positions and type.

    Be Molecule :math:`\\boldsymbol{D}\\in\\mathbb{R}^{(n+1)\\times m}` with
    spatial dimensions :math:`i=1,\\dots,n` and atom number
    :math:`j=1,\\dots,m`. The data matrix then has following structure

    .. math::

        \\boldsymbol{D}=\\begin{bmatrix}
        \\begin{pmatrix}d_{11}\\\\d_{12}\\\\\\vdots\\\\d_{1m}\\end{pmatrix}&
        \\begin{pmatrix}d_{21}\\\\d_{22}\\\\\\vdots\\\\d_{2m}\\end{pmatrix}&
        \\dots&
        \\begin{pmatrix}d_{n1}\\\\d_{n2}\\\\\\vdots\\\\d_{nm}\\end{pmatrix}&
        \\begin{pmatrix}d_{t1}\\\\d_{t2}\\\\\\vdots\\\\d_{tm}\\end{pmatrix}
        \\end{bmatrix}

    with atom data :math:`d_{ij}` of atom :math:`j` and chemical atom type entry
    :math:`t:=n+1`.

    The requested dimension is defined in the parameter module.
    Most of the following functions have been written for :math:`n`-dimensions.
    However, they only have been tested for the three-dimensional case.

    Functions have been provided for editing, moving and transforming the
    entries of the data matrix either for the whole matrix or for specified
    atoms. Private functions are solely needed for internal purposes,
    whereas the public functions are intended to be used for the editing.

    Parameters
    ----------
    name : string, optional
        Molecule name
    short : string, optional
        Molecule short name
    inp : None, string, list, optional
        None for empty Molecule, string to read a molecule from a
        specified file link or a list of either molecules to concatenate these
        into one object, or a data list in the earlier discussed format

    Examples
    --------
    Following example generates a benzene molecule without hydrogen atoms

    .. code-block:: python

        from porems.molecule import Molecule

        mol = Molecule("benzene", "BEN")
        mol.add("C", [0,0,0])
        mol.add("C", 0, r=0.1375, theta= 60)
        mol.add("C", 1, r=0.1375, theta=120)
        mol.add("C", 2, r=0.1375, theta=180)
        mol.add("C", 3, r=0.1375, theta=240)
        mol.add("C", 4, r=0.1375, theta=300)
    """
    def __init__(self, name="Molecule", short="MOL", inp=None):
        # Initialize
        self._box = None
        self._masses = None
        self._mass = None
        self._com = None
        self._mol_list = [self]

        self._dim = 3
        self._charge = 0

        self._name = name
        self._short = short

        # Check data input
        if inp is None:
            self._data = [[] for i in range(self._dim+1)]
        else:
            # Read from file
            if isinstance(inp, str):
                self._data = self._read(inp, inp.split(".")[-1].upper())
            # Concat multiple molecules
            elif isinstance(inp, list):
                # Data list is provided
                if(isinstance(inp[0], list)):
                    self._data = inp
                # List of molecules is provided
                else:
                    self._data = self._concat(inp)


    ################################
    # Private Methods - Management #
    ################################
    def _read(self, file_path, file_type):
        """Read a molecule from a file. Currently only **GRO**, **PFB** and
        **MOL2** files are supported. In case the ``file_type`` is **OBJ**, the
        function will attempt unpickling the object file. Note that loading the
        object this way, only the data will be imported.

        Parameters
        ----------
        file_path : string
            Link to requested file
        file_type : string
            Int for types list id or str for file extension name

        Returns
        -------
        data : list
            data matrix
        """
        # Process input
        if not file_type in ["GRO", "PDB", "MOL2", "OBJ"]:
            print("Unsupported filetype.")
            return

        # Load object
        if file_type == "OBJ":
            try:
                return utils.load(file_path).get_data()
            except:
                print("Erroneous object file.")
                return

        # Read molecule
        data = []
        with open(file_path, "r") as file_in:
            line_idx = 0
            for line in file_in:
                line_val = line.split()

                # Gro file
                if file_type == "GRO":
                    if line_idx > 0 and len(line_val) > 3:
                        coord = [float(line_val[i]) for i in range(3, 5+1)]
                        coord.append(''.join([i for i in line_val[1] if not i.isdigit()]))
                        data.append(coord)

                # Pdb file
                elif file_type == "PDB":
                    if line_val[0] in ["ATOM", "HETATM"]:
                        coord = [float(line_val[i])/10 for i in range(6, 8+1)]
                        coord.append(line_val[11])
                        data.append(coord)

                # Mol2 file
                elif file_type == "MOL2":
                    if len(line_val) > 8:
                        coord = [float(line_val[i])/10 for i in range(2, 4+1)]
                        coord.append(''.join([i for i in line_val[1] if not i.isdigit()]))
                        data.append(coord)

                # Running variable
                line_idx += 1

        # Transform to column
        return utils.column(data)

    def _concat(self, mol_list):
        """Concatenate a molecule list into one molecule object.

        Parameters
        ----------
        mol_list : list
            List of molecule objects to be concatenated

        Returns
        -------
        data : list
            data matrix
        """
        # Initialize
        data = [[] for x in range(self._dim+1)]

        # Concatenate
        for mol in mol_list:
            # Molecules
            mol_data = mol.get_data()
            for i in range(len(mol_data)):
                for j in range(len(mol_data[i])):
                    data[i].append(mol_data[i][j])

        return data

    def _temp(self, atoms):
        """Create a temporary molecule of specified atom ids.

        Parameters
        ----------
        atoms : list
            List of atoms to be duplicated

        Returns
        -------
        mol : Molecule
            Molecule object
        """
        # Create temp matrix
        temp = [[] for x in range(self._dim+1)]

        # Fill matrix
        for i in range(len(atoms)):
            for j in range(self._dim+1):
                temp[j].append(self._data[j][atoms[i]])

        # Create molecule
        return Molecule(inp=temp)

    def _update(self, data, atoms):
        """Exchange data of specified atom with given data.

        Parameters
        ----------
        data : list
            New data list
        atoms : list
            List of atoms to be updated
        """
        for i in range(len(atoms)):
            for j in range(self._dim):
                self._data[j][atoms[i]] = data[j][i]

    def _append(self, mol):
        """Append a given molecule to the current object.

        Parameters
        ----------
        mol : Molecule
            Molecule object
        """
        temp = mol.get_data()
        for i in range(len(temp)):
            for j in range(len(temp[i])):
                self._data[i].append(temp[i][j])


    ##############################
    # Private Methods - Geometry #
    ##############################
    def _dotproduct(self, vec_a, vec_b):
        """Calculate the dot product of two vectors
        :math:`\\boldsymbol{a},\\boldsymbol{b}\\in\\mathbb{R}^n`

        .. math::

            \\text{dot}(\\boldsymbol{a},\\boldsymbol{b})=
            \\begin{pmatrix}a_1\\\\\\vdots\\\\a_n\\end{pmatrix}\\cdot
            \\begin{pmatrix}b_1\\\\\\vdots\\\\b_n\\end{pmatrix}=
            a_1\\cdot b_1+a_2\\cdot b_2+\\dots+a_n\\cdot b_n.

        Parameters
        ----------
        vec_a : list
            First vector :math:`\\boldsymbol{a}`
        vec_b : list
            Second vector :math:`\\boldsymbol{b}`

        Returns
        -------
        dot : float
            Dot product value
        """
        return sum((a*b) for a, b in zip(vec_a, vec_b))

    def _length(self, vec):
        """Calculate the length of a vector
        :math:`\\boldsymbol{a}\\in\\mathbb{R}^n`

        .. math::

            \\text{length}(\\boldsymbol{a})=|\\boldsymbol{a}|
            =\\sqrt{\\boldsymbol{a}\cdot\\boldsymbol{a}}

        Parameters
        ----------
        vec : list
            Vector a

        Returns
        -------
        length : float
            Vector length
        """
        return math.sqrt(self._dotproduct(vec, vec))

    def _vector(self, pos_a, pos_b):
        """Calculate the vector between to two positions
        :math:`\\boldsymbol{a},\\boldsymbol{b}\\in\\mathbb{R}^n`

        .. math::

            \\text{vec}(\\boldsymbol{a},\\boldsymbol{b})
            =\\begin{pmatrix}b_1-a_1\\\\\\vdots\\\\b_n-a_n\\end{pmatrix}

        The two inputs can either be atom indices or to vectoral positions.

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
        elif not isinstance(pos_a, list) and isinstance(pos_b, list):
            print("Vector: Wrong input...")
            return

        # Check dimensions
        if not len(pos_a) == len(pos_b) or not len(pos_a) == self._dim:
            print("Vector: Wrong dimensions...")
            return

        # Calculate vector
        return [pos_b[i]-pos_a[i] for i in range(self._dim)]

    def _unit(self, vec):
        """Transform a vector :math:`\\boldsymbol{a}\\in\\mathbb{R}^n` into a
        unit vector

        .. math::

            \\text{unit}(\\boldsymbol{a})
            =\\frac{\\boldsymbol{a}}{|\\boldsymbol{a}|}

        Parameters
        ----------
        vec : list
            Vector a

        Returns
        -------
        vec : list
            Vector
        """
        length = self._length(vec)

        return [vec[i]/length if not length == 0 else vec[i] for i in range(self._dim)]

    def _cross(self, vec_a, vec_b):
        """Calculate the cross product of two three-dimensional vectors
        :math:`\\boldsymbol{a},\\boldsymbol{b}\\in\\mathbb{R}^3`

        .. math::

            \\text{cross}(\\boldsymbol{a},\\boldsymbol{b})=\\begin{pmatrix}
            a_2\\cdot b_3-a_3\\cdot b_2\\\\
            a_3\\cdot b_1-a_1\\cdot b_4\\\\
            a_1\\cdot b_2-a_2\\cdot b_1
            \\end{pmatrix}

        Parameters
        ----------
        vec_a : list
            First vector :math:`\\boldsymbol{a}`
        vec_b : list
            Second vector :math:`\\boldsymbol{b}`

        Returns
        -------
        vec : list
            Cross product vector
        """
        vec = []
        vec.append(vec_a[1]*vec_b[2]-vec_a[2]*vec_b[1])
        vec.append(vec_a[2]*vec_b[0]-vec_a[0]*vec_b[2])
        vec.append(vec_a[0]*vec_b[1]-vec_a[1]*vec_b[0])

        return vec

    def _angle_polar(self, pos, is_deg=False):
        """Calculate the polar angle of a position vector
        :math:`\\boldsymbol{a}\\in\\mathbb{R}^3`, which is the angle of the
        x-axis towards the reflected position vector on the x-y-plane

        .. math::

            \\text{polar}(\\boldsymbol{a})=\\arctan2(x,y)\\left\\{
            \\begin{array}{ll}
            \\tan^{-1}\\left(\\frac{y}{x}\\right)&x>0\\\\
            \\tan^{-1}\\left(\\frac{y}{x}\\right)+\\pi&x<0,y>0\\\\
            \\pm\\pi&x<0,y=0\\\\
            \\tan^{-1}\\left(\\frac{y}{x}\\right)-\\pi&x<0,y<0\\\\
            +\\frac{\\pi}{2}&x=0,y>0\\\\
            -\\frac{\\pi}{2}&x=0,y<0
            \\end{array}
            \\right.

        with :math:`x` as the first vector entry and :math:`y` as the second.

        Parameters
        ----------
        pos : list
            Position vector :math:`\\boldsymbol{a}`
        is_deg : bool, optional
            True if the output should be in degree

        Returns
        -------
        angle : float
            Polar angle
        """
        try:
            angle = math.atan2(pos[1], pos[0])
        except:
            angle = math.atan(0)

        return angle*180/math.pi if is_deg else angle

    def _angle_azi(self, pos, is_deg=False):
        """Calculate the azimuthal angle of a position vector
        :math:`\\boldsymbol{a}\\in\\mathbb{R}^3`, which is the angle of the
        position vector towards the x-y-plane

        .. math::

            \\text{azimut}(\\boldsymbol{a})
            =\\cos^{-1}\\frac{y}{|\\boldsymbol{a}|}

        with :math:`y` as the second vector entry.

        Parameters
        ----------
        pos : list
            Position vector :math:`\\boldsymbol{a}`
        is_deg : bool, optional
            True if the output should be in degree

        Returns
        -------
        angle : float
            Azimuthal angle
        """
        try:
            angle = math.acos(pos[2]/self._length(pos))
        except:
            angle = math.acos(0)

        return angle*180/math.pi if is_deg else angle

    def _axis(self, inp):
        """Return the three-dimensional unit-vector of the main axes.
        Input is either integer or string

        * 1 or "x" for the x-axis
        * 2 or "y" for the y-axis
        * 3 or "z" for the z-axis

        Parameters
        ----------
        inp : integer, string
            Axis type input

        Returns
        -------
        vec : list
            Unit vector
        """
        # Error message
        axis_error = "Wrong axis definition..."

        # Process input
        if isinstance(inp, str):
            if inp == "x":
                axis = 1
            elif inp == "y":
                axis = 2
            elif inp == "z":
                axis = 3
            else:
                return axis_error
        elif isinstance(inp, int):
            if inp == 1 or inp == 2 or inp == 3:
                axis = inp
            else:
                return axis_error
        else:
            return axis_error

        # Return vector
        return [1 if i == axis-1 else 0 for i in range(self._dim)]

    def _focal(self):
        """Calculate the positional focal-point

        .. math::

            \\text{focal}=\\begin{pmatrix}f_1&f_2&\\dots&f_n\\end{pmatrix}^T

        with

        .. math::

            f_i=\\frac{1}{m}\\sum_j^m d_{ij}

        Returns
        -------
        focal : list
            Focal point
        """
        return [sum(self._data[i])/len(self._data[i]) for i in range(self._dim)]

    def _rotate(self, data, axis, angle, is_deg):
        """Rotate a vector :math:`\\boldsymbol{a}\\in\\mathbb{R}^3`
        along an axis :math:`\\boldsymbol{b}\\in\\mathbb{R}^3` with angle
        :math:`\\alpha\\in\\mathbb{R}`.
        Input for the axis is either a vector or the input for
        function :func:`_vector`.
        The rotation is performed using the rotation-matrix

        .. math::

            \\boldsymbol{R}_\\boldsymbol{n}(\\alpha)=\\begin{pmatrix}
            n_1^2 (1-\\cos\\alpha)+   \\cos\\alpha&n_1n_2(1-\\cos\\alpha)-n_3\\sin\\alpha&n_1n_3(1-\\cos\\alpha)+n_2\\sin\\alpha\\\\
            n_2n_1(1-\\cos\\alpha)+n_3\\sin\\alpha&n_2^2 (1-\\cos\\alpha)+   \\cos\\alpha&n_2n_3(1-\\cos\\alpha)-n_1\\sin\\alpha\\\\
            n_3n_1(1-\\cos\\alpha)-n_2\\sin\\alpha&n_3n_2(1-\\cos\\alpha)+n_1\\sin\\alpha&n_3^2 (1-\\cos\\alpha)+   \\cos\\alpha
            \\end{pmatrix}


        where :math:`n_i` are the entries for the unit vector
        :math:`\\boldsymbol{n}` of the axis. The new coordinates
        :math:`\\boldsymbol{c}` are then calculated using a matrix-vector
        multiplication

        .. math::

            \\boldsymbol{c}=\\boldsymbol{R}_\\boldsymbol{n}\\boldsymbol{a}.

        Parameters
        ----------
        data : list
            Vector :math:`\\boldsymbol{a}`
        axis : integer, string, list
            Axis :math:`\\boldsymbol{b}`
        angle : float
            Angle
        is_deg : bool
            True if the input is in degree

        Returns
        -------
        coord : list
            Vector c as the result of the rotation
        """
        # Angle
        angle = angle*math.pi/180 if is_deg else angle

        # Set vector
        if isinstance(axis, list):
            if len(axis) == self._dim:
                n = axis
            elif len(axis) == 2:
                n = self._vector(axis[0], axis[1])
            else:
                print("Wrong vector dimensions.")
                return
        else:
            n = self._axis(axis)
            if isinstance(n, str):
                print(n)
                return

        # Unit vector
        n = self._unit(n)

        # Define rotation matrix
        n1 = n[0]
        n2 = n[1]
        n3 = n[2]

        c = math.cos(angle)
        s = math.sin(angle)

        r = [[n1*n1*(1.-c) + c, n1*n2*(1.-c) - n3*s,  n1*n3*(1.-c) + n2*s],
             [n2*n1*(1.-c) + n3*s, n2*n2*(1.-c) + c,  n2*n3*(1.-c) - n1*s],
             [n3*n1*(1.-c) - n2*s, n3*n2*(1.-c) + n1*s,  n3*n3*(1.-c) + c]]

        # Rotate
        return [data[0]*r[i][0]+data[1]*r[i][1]+data[2]*r[i][2] for i in range(self._dim)]

    def _box_size(self):
        """Calculate the box size of the current molecule. This is done by
        determining the maximal coordinate value of all atoms in all dimensions

        .. math::

            \\boldsymbol{b}=\\begin{pmatrix}\\max(\\boldsymbol{d}_1)&max(\\boldsymbol{d}_1)&\\dots&max(\\boldsymbol{d}_n)\\end{pmatrix}^T

        where :math:`\\boldsymbol{d}_i` is the dimension-vector of the data
        matrix.

        Returns
        -------
        box : list
            Box length of the current molecule
        """
        return [max(self._data[i]) if max(self._data[i]) > 0 else 0.001 for i in range(self._dim)]


    ##############################
    # Public Methods - Transform #
    ##############################
    def translate(self, vec):
        """Translate data matrix :math:`\\boldsymbol{D}` along a vector
        :math:`\\boldsymbol{a}\\in\\mathbb{R}^n`.

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
        # Translate coordinates
        for i in range(self._dim):
            for j in range(len(self._data[i])):
                self._data[i][j] += vec[i]

    def rotate(self, axis, angle, is_deg=True):
        """Rotate data matrix :math:`\\boldsymbol{D}` around an axis
        :math:`\\boldsymbol{a}\\in\\mathbb{R}^3` with angle
        :math:`\\alpha\\in\\mathbb{R}` using rotation function :func:`_rotate`.

        Parameters
        ----------
        axis : integer, string, list
            Axis
        angle : float
            Angle
        is_deg : bool, optional
            True if the input is in degree
        """
        # Initialize
        dim = self._dim
        data = self._data

        # Rotate
        for i in range(len(data[0])):
            coord = [0 for i in range(dim)]
            for j in range(dim):
                coord[j] += data[j][i]

            coord = self._rotate(coord, axis, angle, is_deg)
            for j in range(dim):
                data[j][i] = coord[j]

    def put(self, atom, pos):
        """Change the position of an atom to a given vector
        :math:`\\boldsymbol{a}\\in\\mathbb{R}^n`.

        Parameters
        ----------
        atom : integer, string
            Atom-id whose position will be changed, use **"last"** for the
            newest atom
        pos : list
            New position vector
        """
        # User input
        atom = len(self._data[0])-1 if atom == "last" else atom

        # Change position
        for i in range(len(pos)):
            self._data[i][atom] = pos[i]

    def move(self, atom, pos):
        """Move whole the molecule to a new position, where the dragging point
        is a given atom that is moved to a specified position.

        Parameters
        ----------
        atom : integer, string
            Main atom-id whose position will be changed, use **"last"** for the
            newest atom
        pos : list
            New position vector
        """
        # Calculate vector
        if atom == "last":
            atom = len(self._data[0])-1

        vec = self._vector(self.pos(atom), pos)

        # Translate
        self.translate(vec)

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
        vec = [pos[i]-min(self._data[i]) for i in range(self._dim)]

        # Resett box size
        self._box = None

        # Translate molecule
        self.translate(vec)

        return vec

    def box(self, size, is_extend=True):
        """Create a molecule box around a centred molecule, where the
        coordinate systems zero point is set to the positional focal point of
        the molecule. The resulting box has the edge length of the given size
        value. If ``is_extend`` is set to True, then the maximal molecule
        coordinate of all dimensions is added to ``size``.

        Parameters
        ----------
        size : float
            Box size
        is_extend : bool, optional
            Should the molecule size be added to the box-size
        """
        # Get focal point
        focal = self._focal()

        # Calculate box size
        boxFocal = [f+size for f in focal]
        self._box = [max(boxFocal) if is_extend else size for f in range(3)]

        # Move molecule to box center
        self.add("R", focal)
        self.move("last", [b/2 for b in self._box])
        self.delete("last")


    #########################
    # Public Methods - Edit #
    #########################
    def part_move(self, bond, atoms, length, vec=None):
        """Change the length of a specified bond. Variable ``atoms`` specifies
        which atoms or rather which part of the molecule needs to be moved for
        this specific bond. The given length is going to be the new bond length,
        not by how much the bond length is changed.

        The move vector is determined automatically by the given length and atom
        bond. This vector can also be given manually with no regards to length,
        by setting the vector ``vec``.

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
            Set this vector, to manually move the atoms

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
        length = abs(length-self.bond(*bond)[1])

        # Set vector
        if vec == None:
            vec = self._vector(bond[0], bond[1])
        vec = [v*length for v in self._unit(vec)]

        # Move molecule
        temp.translate(vec)

        # Update positions
        self._update(temp.get_data(), atoms)

    def part_rotate(self, bond, atoms, angle, zero):
        """Rotate a set of specified atoms around a given bond as the rotation
        axis. First however the system needs to be set to zero. Therefore the
        atom id to define the new coordinate system must be given for the set
        of specified atoms. Normally this is the atoms that lies at on an end
        of the given bond axis.

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

        # Update positions
        self._update(temp.get_data(), atoms)

    def part_angle(self, bond_a, bond_b, atoms, angle, zero):
        """Change the bond angle of two bond vectors. Variable ``atoms``
        specifies which atoms or rather which part of the molecule needs to be
        rotated in order to change the specified bond angle. First however the
        system needs to be set to zero. Therefore, the atom id to define the new
        coordinate system has to be given for the set of specified atoms.
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
        if len(bond_a) == 2 and len(bond_b) == 2:
            vec = self._cross(self._vector(bond_a[0], bond_a[1]),
                              self._vector(bond_b[0], bond_b[1]))
        elif len(bond_a) == self._dim and len(bond_b) == self._dim:
            vec = self._cross(bond_a, bond_b)
        else:
            print("Wrong bond input...")
            return
        temp.rotate(vec, angle)

        # Update positions
        self._update(temp.get_data(), atoms)

    def length_angle(self, bond, bond_a, bond_b, atoms, zero, length, angle=[0, 45], grid=0.001, is_silent=True, is_negative=False):
        """If the bond length of an exemplary circular molecule is correlated
        to an angle, then the bond length has to be changed in dependence to
        the bond angle.

        Variable ``atoms`` specifies which atoms or rather which part of the
        molecule needs moved and rotated to change the specified bond length.

        First however the system needs to be set to zero. Therefore, the atom id
        to define the new coordinate system must be given for the set of
        specified atoms. Normally this is the atom that touches the angle.

        The rotation axis is determined by creating the cross product
        of the two bond vectors. Thus, getting the normal vector of a surface
        that contains both bond vectors.

        The new length is determined, by rotating the specified molecule part in
        steps and calculating the resulting length error. If the error is zero,
        then the rotation steps stop.

        Parameters
        ----------
        bond : list
            Bond of which the length should be changed
        bond_a : list
            First bond vector given as a list of two atom ids
        bond_b : list
            Second bond vector given as a list of two atom ids
        atoms : integer, list
            List of atoms to be rotated (can also be one id)
        zero : integer
            Atom id to define zero point of the new coordinate system
        length : float
            New bond length
        angle : list, optional
            Angle range to be tested
        grid : float, optional
            angle steps fineness
        is_silent : bool, optional
            True to suppress error-value messages
        is_negative : bool, optional
            True for a negative rotation

        Examples
        --------
        .. code-block:: python

            mol.length_angle([0, 19], [8, 9], [9, 10], [10, 11, 12], mol.pos(9))

        """
        angles = [x*grid for x in range(int(angle[0]/grid), int(angle[1]/grid)+1)]
        for a in angles:
            a = -a if is_negative else a
            self.part_angle(bond_a, bond_b, atoms, a, zero)

            error = round(abs(self.bond(*bond)[1]-length), 3)

            if not is_silent:
                print("Error = "+str(error))

            if error == 0:
                break


    ##############################
    # Public Methods - Calculate #
    ##############################
    def angle(self, vec_a, vec_b, is_deg=True):
        """Calculate the angle between two vectors
        :math:`\\boldsymbol{a},\\boldsymbol{b}\\in\\mathbb{R}^n`

        .. math::

            \\text{angle}=\\cos^{-1}\\frac{\\boldsymbol{a}\cdot\\boldsymbol{b}}
            {|\\boldsymbol{a}||\\boldsymbol{a}|}

        Parameters
        ----------
        vec_a : list
            First vector :math:`\\boldsymbol{a}`
        vec_b : list
            Second vector :math:`\\boldsymbol{b}`
        is_deg : bool, optional
            True if the output should be in degree

        Returns
        -------
        angle : float
            Angle
        """
        angle = math.acos(self._dotproduct(vec_a, vec_b)/(self._length(vec_a)*self._length(vec_b)))

        return angle*180/math.pi if is_deg else angle

    def pos(self, atom):
        """Get the position of an atom as a vector.

        Parameters
        ----------
        atom : integer
            Atom id

        Returns
        -------
        pos : list
            Position vector of the specified atom
        """
        return [self._data[i][atom] for i in range(self._dim)]

    def bond(self, inp_a, inp_b):
        """Return the bond length and vector of a specified bond. The two inputs
        can either be atom indices or to vectoral positions.

        Parameters
        ----------
        inp_a : integer, list
            Either an atom id or a position vector
        inp_b : integer, list
            Either an atom id or a position vector

        Returns
        -------
        bond : list
            Bond vector and length

        Examples
        --------
        .. code-block:: python

            mol.bond(0, 1)
            mol.bond(*[0, 1])
            mol.bond([1, 0, 0], [0, 0, 0])
        """
        return [self._vector(inp_a, inp_b), self._length(self._vector(inp_a, inp_b))]


    #########################
    # Public Methods - Edit #
    #########################
    def add(self, name, pos, bond=None, r=0, theta=0, phi=0, is_deg=True):
        """Add a new atom in polar coordinates. The ``pos`` input is either
        an atom id, that determines is the bond-start,
        or a vector for a specific position.

        Bond has to be given, if the polar coordinates are dependent on the bond
        vector as the basic axis. The coordinate system is then transformed
        to the new bond axis. If set to None, then the given coordinates are
        thought to be dependent on the basic axes.

        If the given position is an atom id, then a bond is automatically added
        to the bond matrix.

        Parameters
        ----------
        name : string
            Atom type
        pos : integer, list
            Position of the atom
        bond : None, list, optional
            Bond axis
        r : float, optional
            Bond length
        theta : float, optional
            Azimuthal angle
        phi : float, optional
            Polar angle
        is_deg : bool, optional
            True if the input of the angles in degree

        Examples
        --------
        .. code-block:: python

            mol.add("C", [0, 0, 0])
            mol.add("C", 0, r=0.153, theta=-135)
            mol.add("C", 1, [0, 1], r=0.153, theta= 135)
        """
        # Angles
        phi *= math.pi/180 if is_deg else 1
        theta *= math.pi/180 if is_deg else 1

        # Transform spherical to cartesian coordinates
        x = r*math.sin(theta)*math.cos(phi)
        y = r*math.sin(theta)*math.sin(phi)
        z = r*math.cos(theta)
        coord = [x, y, z]

        # Bond vector
        if bond == None:
            vec = self._axis("z")
        else:
            vec = self._vector(bond[0], bond[1])

        # Calculate angles for rotation
        phi = self._angle_polar(vec)
        theta = self._angle_azi(vec)

        # Rotate towards axis
        tRot = self._cross([0, 0, 1], vec)
        if sum(tRot) == 0:
            tRot = "y"

        coord = self._rotate(coord, tRot, theta, is_deg=False)

        # Process position input
        b = None
        if isinstance(pos, int):
            b = pos
            pos = self.pos(pos)

        # Add coordinates
        for i in range(self._dim):
            self._data[i].append(pos[i]+coord[i])

        # Add name
        self._data[self._dim].append(name)

    # Delete an atom
    def delete(self, atom):
        """Delete specified atom from the molecule.
        The input can also be a list of atom ids.

        Parameters
        ----------
        atom :  int, list
            Atom id or list to be deleted
        """
        # Process input
        if isinstance(atom, str):
            if atom == "last":
                atom = [len(self._data[0])-1]
        elif not isinstance(atom, list):
            atom = [atom]
        atom = sorted(atom, reverse=True)

        # Remove line
        for data in self._data:
            for a in atom:
                data.pop(a)

    def overlap(self, error=0.005, is_print=False):
        """Check if atoms are overlapping (as in duplicates) and delete these.

        Parameters
        ----------
        error : float, optional
            Error of the overlap search
        is_print : bool, optional
            True to print the list of overlapping atoms
        """
        # Initialize
        data = self._data
        dim = self._dim

        # Find duplicate atoms in reverse
        for i in list(reversed(range(self.get_num()))):
            # Get positions
            atom_a = self.pos(i)
            for j in range(self.get_num()):
                if not j == i:
                    atom_b = self.pos(j)

                    # Check if identical
                    sum_del = 0
                    a = []
                    for k in range(dim):
                        if abs(atom_a[k]-atom_b[k]) < error:
                            a.append(abs(atom_a[k]-atom_b[k]))
                            sum_del += 1

                    # Delete atom
                    if sum_del == dim:
                        # print(i,j,a)
                        self.delete(i)
                        if is_print:
                            print(i)
                        break

    def change_type(self, atom, atom_type):
        """Change the atom type of a specified atom.

        Parameters
        ----------
        atom : integer
            Atom id
        atom_type : string
            New atomtype
        """
        self._data[self._dim][atom] = atom_type

    def switch_order(self, atom_a, atom_b):
        """Switch atom order of two atoms.

        Parameters
        ----------
        atom_a : integer
            Atom id of the first atom
        atom_b : integer
            Atom id of the second atom
        """
        # Get data
        data_a = [self._data[i][atom_a] for i in range(self._dim+1)]
        data_b = [self._data[i][atom_b] for i in range(self._dim+1)]

        # Switch data
        for i in range(self._dim+1):
            self._data[i][atom_a] = data_b[i]
            self._data[i][atom_b] = data_a[i]


    ###########################
    # Public Methods - Output #
    ###########################
    def table(self):
        """Create a pandas table of the molecule data."""
        import pandas as pd

        # Create dictionary
        data = {"Axis "+str(dim+1): self._data[dim] for dim in range(self._dim)}
        data["Type"] = self._data[-1]

        # Create dataframe
        return pd.DataFrame(data)


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
        """Set the molecule short name

        Parameters
        ----------
        short : string
            Molecule short name
        """
        self._short = short

    def set_mol_list(self, mol_list):
        """Set the molecule list for writing the structure file

        Parameters
        ----------
        mol_list : list
            List of molecule objects
        """
        self._mol_list = mol_list

    def set_box(self, box):
        """Set the molecule box dimensions.

        Parameters
        ----------
        box : list
            Box dimension
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

    def set_masses(self, masses=None):
        """Set the molar masses of the atoms.

        Parameters
        ----------
        masses : list, None, optional
            List of molar masses in :math:`\\frac g{mol}`
        """
        import porems.database as db
        self._masses = [db.get_mass(atom) for atom in self._data[self._dim]] if masses is None else masses

    def set_mass(self, mass=None):
        """Set the molar mass of the molecule.

        Parameters
        ----------
        mass : float, None, optional
            Molar mass in :math:`\\frac g{mol}`
        """
        self._mass = sum(self.get_masses()) if mass is None else mass

    def set_com(self, com=None):
        """Set the center of mass coordinates of the molecule.

        Parameters
        ----------
        com : list, None, optional
            Centre of mass coordinates
        """
        if com is not None:
            self._com = com
        else:
            self._com = []
            masses = self.get_masses()
            for i in range(self._dim):
                self._com.append(0)
                for j in range(len(self._data[i])):
                    self._com[i] += masses[j]*self._data[i][j]

            self._com = [x/sum(masses) for x in self._com]


    ##################
    # Getter Methods #
    ##################
    def get_data(self):
        """Return the data matrix.

        Returns
        -------
        data : list
            Data matrix
        """
        return self._data

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

    def get_mol_list(self):
        """Return the global molecule list for writing the structure file.

        Returns
        -------
        mol_list : list
            List of all molecules
        """
        return self._mol_list

    def get_box(self):
        """Calculate and return the box size of the current molecule.

        Returns
        -------
        box : list
            Box dimensions
        """
        return self._box_size()

    def get_box_c(self):
        """Return the current box size.

        Returns
        -------
        box : list
            Box dimensions
        """
        return self._box

    def get_num(self):
        """Return the number of atoms.

        Returns
        -------
        num : integer
            Number of atoms
        """
        return len(self._data[0])

    def get_type(self, atom):
        """Return the atomtype of the given atom id.

        Parameters
        ----------
        atom : integer
            Atom id

        Returns
        -------
        box : list
            Box dimensions
        """
        return self._data[self._dim][atom]

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
        if self._masses is None:
            self.set_masses()
        return self._masses

    def get_mass(self):
        """Return the molar mass of the molecule.

        Returns
        -------
        mass : float
            Molar mass in :math:`\\frac g{mol}`
        """
        if self._mass is None:
            self.set_mass()
        return self._mass

    def get_com(self):
        """Return the centre of mass of the molecule.

        Returns
        -------
        com : list
            centre of mass coordinates
        """
        if self._com is None:
            self.set_com()
        return self._com
