################################################################################
# Verlet Class                                                                 #
#                                                                              #
"""Verlet lists definitions."""
################################################################################


import math
import multiprocessing as mp

import porems.utils as utils


class Verlet:
    """This class converts the molecule box into verlet lists.

    The aim of verlet lists is reducing the workload on search algorithms for
    atom pairs. Normally the effort would be

    .. math::

        \\mathcal O(n^2)

    with the number of atoms :math:`n`, since each atom has to be compared with
    each atom if their bond length is within the search distance.
    Since this distance is however confined, the idea is reducing
    the search space by dividing the molecule box into small verlet boxes. Thus
    the atom pairs must only be checked in the box containing the starting atom
    and the small boxes surrounding the main box. In the three-dimensional case,
    the effort reduces to

    .. math::

        \\begin{array}{cc}
            \\mathcal O(n\\cdot N),&N=\\sum_{j=0}^ma_j
        \\end{array}

    with the number of boxes m, in the three dimensional case this is
    :math:`3^3=27`, the number of atoms :math:`a_j` of box :math:`j` and the
    number of surround atoms :math:`N`. This is considerably less and of course
    depending on the verlet box size. Note that the verlet size :math:`v` has to
    be greater than the specified bond length :math:`l` for the search

    .. math::

        v > l.

    Since molecular simulation considers periodic
    boundary conditions, these must also be applied for the surround box search.

    Parameters
    ----------
    mol : Molecule
        Molecule to be divided
    size : float
        Verlet box edge size
    is_pbc : bool
        True if periodic boundary conditions are needed
    """
    def __init__(self, mol, size, is_pbc):
        # Initialize
        self._dim = 3
        self._np = mp.cpu_count()
        self._is_pbc = is_pbc
        self._mol = mol
        self._data = mol.get_data()
        self._box = [[], []]

        self._verlet(size)
        self._fill()


    ################################
    # Private Methods - Definition #
    ################################
    def _verlet(self, size):
        """Here the possible number of verlet boxes is calculated for each
        dimension for the given size and molecule dimensions. The given size is
        adjusted, since only whole numbered box number are reasonable.
        A list of boxes is generated containing the coordinates of the zero
        point of each box. Furthermore, an empty list for each box added, that
        will contain pointer to the atoms.

        The boxes are numbered as follows for the example of a :math:`3\\times3`
        block

        .. math::

            \\begin{array}{ccc}
                \\begin{bmatrix}
                     6& 7& 8\\\\
                     3& 4& 5\\\\
                     0& 1& 2
                \\end{bmatrix}\\text{ for }z=0,&
                \\begin{bmatrix}
                    15&16&17\\\\
                    12&13&14\\\\
                     9&10&11
                \\end{bmatrix}\\text{ for }z=1,\dots
            \\end{array}

        where the x-axis is the first count direction, y-axis the second and
        z-axis the third. The position list for a box id contains the values for
        the x, y and z values. For example the position of box id :math:`12`
        would be

        .. math::

            \\text{pos}(12)=\\begin{bmatrix}0&1&1\\end{bmatrix}.

        Parameters
        ----------
        size : float
            Desired verlet box length
        """
        # Initialize
        boxes = []
        box_size = self._mol.get_box()

        # Calculate box number and sizes
        count = [math.floor(box/size) for box in box_size]
        size = [box_size[i]/count[i] for i in range(len(count))]

        # Fill verlet list - k=x, j=y, i=z
        for i in range(count[2]):
            for j in range(count[1]):
                for k in range(count[0]):
                    boxes.append([size[0]*k, size[1]*j, size[2]*i])

        # Globalize
        self._box[0] = boxes
        self._size = size
        self._count = count

    def _fill(self):
        """Atoms are so to say filled in the verlet boxes depending on their
        position. This is done by adding atom ids to specific boxes in the box
        array.
        """
        # Initialize
        data = self._data
        data = self._mol.get_data()
        dim = self._dim
        size = self._size
        count = self._count
        pointer = [[] for i in range(count[0]*count[1]*count[2])]

        # Fill boxes with atoms
        for i in range(len(data[0])):
            pos = [0 for i in range(dim)]
            for j in range(dim):
                pos[j] = math.floor(data[j][i]/size[j])
                pos[j] -= 1 if pos[j] == count[j] else 0

            pointer[self._index(pos)].append(i)

        # Globalize
        self._box[1] = pointer


    ##############################
    # Private Methods - Iterator #
    ##############################
    def _index(self, pos):
        """Get the box id from a given box position which is a list of x, y and
        z values.

        Parameters
        ----------
        pos : list
            Box position in x, y and z values

        Returns
        -------
        index : integer
            Box index
        """
        # Initialize
        count = self._count
        index = 0

        # Calculate id
        index += pos[0]
        index += pos[1]*count[0]
        index += pos[2]*count[0]*count[1]

        return int(index)

    def _position(self, index):
        """Get the box position in x, y and z values from a given box id.

        Parameters
        ----------
        index : integer
            Box index

        Returns
        -------
        pos : list
            Box position in x, y and z values
        """
        # Initialize
        dim = self._dim
        count = self._count
        pos = [0 for i in range(self._dim)]

        # Get index
        pos[2] = math.floor(index/(count[0]*count[1]))
        pos[1] = math.floor((index-pos[2]*count[0]*count[1])/count[0])
        pos[0] = index-pos[1]*count[0]-pos[2]*count[0]*count[1]

        return pos

    def _input(self, inp):
        """Check if input is a box id or a position. If it is an id then
        the id is transformed to a position list.

        Parameters
        ----------
        inp : integer, list
            Input

        Returns
        -------
        pos : list
            Box position in x,y and z values
        is_id : bool
            True if the input was an id
        """
        if isinstance(inp, list):
            pos = inp
            is_id = False
        elif isinstance(inp, int):
            pos = self._position(inp)
            is_id = True
        else:
            print("Wrong input!")
            return None

        return pos, is_id

    def _forward(self, dim, inp):
        """Helper function for stepping forward with the iterator. Periodic
        boundary conditions are applied depending on the global pbc variable.

        Parameters
        ----------
        dim : integer
            Stepping dimension
        inp : integer, list
            Position input, can also be an id

        Returns
        -------
        pos : integer, list
            Iterated position in the input format
        """
        # Check if input is None
        if inp is None:
            return

        # Process input
        pos, is_id = self._input(inp)

        # Move
        if pos[dim] == self._count[dim]-1:
            pos[dim] = 0 if self._is_pbc else None
        else:
            pos[dim] += 1

        # Return
        if not None in pos:
            return self._index(pos) if is_id else pos

    def _backward(self, dim, inp):
        """Helper function for stepping backward with the iterator. Periodic
        boundary conditions are applied depending on the global pbc variable.

        Parameters
        ----------
        dim : integer
            Stepping dimension
        inp : integer, list
            Position input, can also be an id

        Returns
        -------
        pos : integer, list
            Iterated position in the input format
        """
        # Check if input is None
        if inp is None:
            return

        # Process input
        pos, is_id = self._input(inp)

        # Move
        if pos[dim] == 0:
            pos[dim] = self._count[dim]-1 if self._is_pbc else None
        else:
            pos[dim] -= 1

        # Return
        if not None in pos:
            return self._index(pos) if is_id else pos

    # Iterators
    def _right(self, inp): return self._forward(0, inp)
    def _left(self, inp): return self._backward(0, inp)
    def _top(self, inp): return self._forward(1, inp)
    def _bot(self, inp): return self._backward(1, inp)
    def _back(self, inp): return self._forward(2, inp)
    def _front(self, inp): return self._backward(2, inp)


    ###########################
    # Parallelization methods #
    ###########################
    def _find_bond(self, box_list, atom_type, distance, error, condition):
        """Search for a bond in the given verlet boxes. This function
        automatically searches for the surrounding boxes and checks in all 27 of
        them for the defined bond (atom types and bond length). If the bond
        length is within a given error, then the atoms will be added to an
        output. The user can also define an additional and stricter condition
        for adding the atoms.

        Parameters
        ----------
        box_list : None, list
            None for all verlet boxes, list for specified boxes
        atom_type : list
            List of two atom type stings
        distance : float
            Bond length
        error : float
            Tolerated deviation of the bond length
        condition : None, function
            None or a function f(Molecule, bond) with a boolean output

        Returns
        -------
        bonds : list
            Bond array containing lists of two atom ids
        """
        # Initialize
        data = self._data
        dim = self._dim
        size = self._size
        count = self._count
        box_len = [count[i]*size[i] for i in range(dim)]
        output = []

        # User input
        if box_list == None:
            box_list = [i for i in range(len(self.get_box()[1]))]

        # Loop through all boxes
        for box in box_list:
            # Get surrounding boxes
            neighbour = self.neighbour(box)
            atoms = []
            for n in neighbour:
                atoms.extend(self._box[1][n])

            # Run through atoms in main boxes
            for atom_a in self._box[1][box]:
                # Check type
                if data[3][atom_a] == atom_type[0]:
                    entry = [atom_a, []]
                    pos_a = self._mol.pos(atom_a)
                    # Search in neighbouring boxes for partner
                    for atom_b in atoms:
                        if data[3][atom_b] == atom_type[1] and not atom_a == atom_b:
                            # Calculate length
                            length = self._mol.bond(atom_a, atom_b)[1]

                            # Periodic boundaries
                            if length > max(size)*3 and self._is_pbc:
                                pos_b = self._mol.pos(atom_b)
                                temp = [pos for pos in pos_b]

                                for i in range(dim):
                                    # Without z-axis
                                    if abs(pos_b[i]-pos_a[i]) > size[i]*3:  # and not i==dim-1:
                                        temp[i] += -box_len[i] if pos_b[i] > pos_a[i] else box_len[i]

                                # New length
                                length = self._mol.bond(self._mol.pos(atom_a), temp)[1]

                            # Add if in range
                            if abs(length-distance) < error:
                                entry[1].append(atom_b)

                    # Add atom
                    is_add = True
                    if condition is not None:
                        is_add = condition(self._mol, entry[0])
                    if is_add:
                        output.append(entry)

        return output

    def find_parallel(self, box_list, atom_type, distance, error, condition=None, is_time=False):
        """Parallelized bond search of function :func:`_find_bond`.

        Parameters
        ----------
        box_list : None, list
            None for all verlet boxes, list for specified boxes
        atom_type : list
            List of two atom type stings
        distance : float
            Bond length
        error : float
            Tolerated deviation of the bond length
        condition : None, function, optional
            None or a function f(Molecule, bond) with a boolean output
        is_time : bool, optional
            True to print the used search time

        Returns
        -------
        bonds : list
            Bond array containing lists of two atom ids
        """
        # Initialize
        np = self._np

        # User input
        if box_list == None:
            box_list = [i for i in range(len(self.get_box()[1]))]

        # Time calculation
        if is_time:
            t = m.tic()

        # Break box list down on number of processors
        box_len = math.floor(len(box_list)/np)

        boxes = []
        for i in range(np):
            if i == np-1:
                boxes.append(box_list[box_len*i:])
            else:
                boxes.append(box_list[box_len*i:box_len*(i+1)])

        # Parallel search
        pool = mp.Pool(processes=np)
        results = [pool.apply_async(self._find_bond, args=(
            x, atom_type, distance, error, condition)) for x in boxes]

        # Concatenate processes
        output = []
        for p in results:
            output.extend(p.get())

        # Time calculation
        if is_time:
            m.toc(t, "Find bond "+atom_type[0]+"-"+atom_type[1])

        return output


    ##################
    # Public methods #
    ##################
    def neighbour(self, index, is_self=True):
        """Get the ids of the boxes surrounding the given one.

        Parameters
        ----------
        index : integer
            Main box id
        is_self : bool, optional
            True to add the main box to the output

        Returns
        -------
        neighbour : list
            List of surrounding boxes
        """
        # Initialize
        neighbour = []

        # Find neighbours
        z = [self._back(index), index, self._front(index)]
        y = [[self._top(i), i, self._bot(i)] for i in z]

        for i in range(len(z)):
            for j in range(len(y[i])):
                neighbour.append(self._left(y[i][j]))
                neighbour.append(y[i][j])
                neighbour.append(self._right(y[i][j]))

        if not is_self:
            neighbour.pop(13)

        return [n for n in neighbour if n is not None]

    def delete(self, inp):
        """Delete given atoms, then refill the verlet boxes,
        since the list ids changed.

        Parameters
        ----------
        inp : list, integer
            Atom list or single atom id
        """
        self._mol.delete(inp)
        self._fill()


    ##################
    # Setter methods #
    ##################
    def set_pbc(self, pbc):
        """Turn the periodic boundary conditions on and off.

        Parameters
        ----------
        bond : bool
            PBC mode
        """
        self._is_pbc = pbc


    ##################
    # Getter methods #
    ##################

    def get_box(self):
        """Return the box list.

        Returns
        -------
        box : list
            Box list
        """
        return self._box

    def get_count(self):
        """Return the number of boxes in each dimension.

        Returns
        -------
        count : list
            Box number
        """
        return self._count

    def get_size(self):
        """Return the box size in each dimension.

        Returns
        -------
        size : list
            Box size
        """
        return self._size

    def get_mol(self):
        """Return the molecule.

        Returns
        -------
        mol : Molecule
            Molecule
        """
        return self._mol
