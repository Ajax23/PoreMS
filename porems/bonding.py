################################################################################
# Bonding Class                                                                #
#                                                                              #
"""Bonding search structure."""
################################################################################


import math

import porems.utils as utils


class Bonding:
    """The aim of this class is preserving all information of the silicon grid
    bonds. This is needed since firstly, verlet searches get costly the more
    often they are needed and secondly, multiple searches deteriorate the bond
    information, due to numeric errors thus resulting in bonds being lost.

    The idea here was reducing the verlet searches to a single search by
    creating a connection matrix of all oxygen and silicon atoms.
    In fact two matrices :math:`ox\\in\\mathbb{N}^{no\\times2}` for oxygen and
    :math:`si\\in\\mathbb{N}^{ns\\times2}` were created with numbers of
    oxygen atoms :math:`no` and silicon atoms :math:`ns`

    .. math::

        \\begin{array}{cc}
            ox = \\left[
                \\begin{array}{ll}
                    o_0&[s_{0,0},s_{0,1}]\\\\
                    o_1&[s_{1,0},s_{1,1}]\\\\
                    &\\vdots\\\\
                    o_{no}&[s_{no,0},s_{no,1}]
                \\end{array}
            \\right],&
            si = \\left[
                \\begin{array}{ll}
                    s_0&[o_{0,0},o_{0,1},o_{0,2},o_{0,3}]\\\\
                    s_1&[o_{1,0},o_{1,1},o_{1,2},o_{1,3}]\\\\
                    &\\vdots\\\\
                    s_{ns}&[o_{ns,0},o_{ns,1},o_{ns,2},o_{ns,3}]
                \\end{array}
            \\right]
        \\end{array}

    with oxygen :math:`o_{i=0,\\dots,no}` and silicon :math:`s_{j=0,\\dots,ns}`.
    The entries next to the main atoms are the silicon atoms
    :math:`s_{i,ko=0,\\dots,so}` bound to oxygen :math:`i` and oxygen atoms
    :math:`o_{j,ks=0,\\dots,os}` bound to silicon :math:`j`, which are found
    using the verlet search algorithm. As an optimization, these values are not
    the atom ids, but the list pointers of the specific atoms in their
    corresponding connection matrices :math:`si` and :math:`ox`. Due to chemical
    properties, the maximal possible silicon atoms bound to one oxygen is
    :math:`so=2` and the maximal possible number of oxygen atoms bound to one
    silicon is :math:`os=4`.

    In the beginning, considering periodic boundary conditions and before
    removing any atoms, every atom should be saturated with partners. This
    means that every silicon has a set of four oxygens and every oxygen has two
    silicon bonds.

    With this connection matrix, bonds can be easily deleted and their number
    for specific atoms can be easily determined without needing another search
    call.

    Parameters
    ----------
    verlet : Verlet
        Verlet list object
    """

    def __init__(self, verlet):
        # Initialize
        osi = []
        sio = []
        pointer = []

        self._mol = verlet.get_mol()
        self._data = self._mol.get_data()
        self._dim = 3
        self._si_grid = 0.155
        self._y_grid = 0.073

        # Search for bounds
        bonds = verlet.find_parallel(None, ["Si", "O"], self._si_grid, 10e-2)

        # Fill atoms in list
        for i in range(self._mol.get_num()):
            atom_type = self._mol.get_type(i)
            if atom_type == "Si":
                pointer.append(len(sio))
                sio.append([i, []])
            elif atom_type == "O":
                pointer.append(len(osi))
                osi.append([i, []])

        # Transform to column
        sio = utils.column(sio)
        osi = utils.column(osi)

        # Add bonds to lists
        for i in range(len(bonds)):
            ids = pointer[bonds[i][0]]
            ido = [pointer[bonds[i][1][j]] for j in range(len(bonds[i][1]))]

            sio[1][ids] = ido
            for o in ido:
                osi[1][o].append(ids)

        # Make global
        self._sio = sio
        self._osi = osi
        self._pointer = pointer
        self._remove = []


    ###################
    # Private Methods #
    ###################
    def _unbind(self, silicon, oxygen):
        """Remove the bond between a silicon and an oxygen atom from the
        connection matrix.

        Parameters
        ----------
        silicon : integer
            Silica atom id
        oxygen : integer
            Oxygen atom id
        """
        # initialize
        sio = self._sio
        osi = self._osi
        pointer = self._pointer

        # Get pointer
        id_s = pointer[silicon]
        id_o = pointer[oxygen]

        # Remove bond
        sio[1][id_s].pop(sio[1][id_s].index(id_o))
        osi[1][id_o].pop(osi[1][id_o].index(id_s))

    def _geminal(self):
        """Get the list silicon that have geminal bonds - two binding sites.

        Returns
        -------
        geminal : list
            List of silicon with geminal bonds
        """
        # Get list of geminal silicon
        oxygen = [o[0] for o in self._osi[1] if len(o) == 1]
        silicon = list(set(oxygen))
        geminal = [self._sio[0][s] for s in silicon if oxygen.count(s) == 2]

        return geminal


    ############################
    # Public Methods - Editing #
    ############################
    def unlink(self, atoms):
        """Remove all bonds of a specified atom from the connection matrix.

        Parameters
        ----------
        atoms : list, integer
            List of atoms to be unlinked, can also be one atom id
        """
        # initialize
        sio = self._sio
        osi = self._osi
        pointer = self._pointer
        data = self._data

        # User input
        if not isinstance(atoms, list):
            atoms = [atoms]

        # Remove atoms
        for atom in atoms:
            # Get type
            atom_type = data[self._dim][atom]
            if atom_type == "Si":
                list_a = sio
                list_b = osi
            elif atom_type == "O":
                list_a = osi
                list_b = sio

            # Get indices
            id_a = pointer[atom]
            atoms_b = list_a[1][id_a]

            # Unlink
            for ab in atoms_b:
                list_b[1][ab].pop(list_b[1][ab].index(id_a))

            list_a[1][id_a] = []

    def remove(self, atoms):
        """Add an atom to a list of atoms to be deleted. This method is
        important since by this call no bonds are removed.

        Parameters
        ----------
        atoms : list, integer
            List of atoms to be removed, can also be one atom id
        """
        # User input
        if not isinstance(atoms, list):
            atoms = [atoms]

        # Add to remove list
        self._remove.extend(atoms)

    def delete(self, is_list=False):
        """Delete all atoms without any bonds and the atoms in the remove list.

        Parameters
        ----------
        is_list : bool, optional
            True to return the atom list to be deleted, False to delete the
            atoms with no return value.

        Returns
        -------
        del : list
            List of atoms to be deleted
        """
        # Initialize
        sio = self._sio
        osi = self._osi

        # Search for empty bonds
        del_list = []
        lists = [sio, osi]

        for target in lists:
            for i in range(len(target[0])):
                if len(target[1][i]) == 0:
                    del_list.append(target[0][i])

        # Concatenate remove list
        del_list.extend(self._remove)
        del_list = list(set(del_list))

        # Return list or delete
        return del_list if is_list else self._mol.delete(del_list)


    ##############################
    # Public Methods - Structure #
    ##############################
    def attach(self):
        """The bonds in drilling direction are broken, since on these sites the
        reservoirs containing molecules are appended and binding sites are
        needed. After breaking the bonds, oxygens are added to the unsaturated
        silicon on the outer surface.
        """
        # Initialize
        sio = self._sio
        osi = self._osi
        pointer = self._pointer
        data = self._data
        mol = self._mol
        si_grid = self._si_grid
        y_grid = self._y_grid

        drill = mol.get_drill()
        gap = mol.get_gap()
        box = mol.get_box_c()

        angle = 2*math.atan(math.sqrt(2))*180/math.pi

        # Atom adding helper function
        def add_oxy(si, o, dist, theta=0):
            # Process input
            if not isinstance(o,    list):
                o = [o]
            if not isinstance(theta, list):
                theta = [theta]

            # Run through oxygens
            for i in range(len(o)):
                # Break bond
                self._unbind(sio[0][si], osi[0][o[i]])

                # Add atom
                mol.add("O", mol.pos(sio[0][si]), r=dist, theta=theta[i])

                # Add bond
                osi[0].append(mol.get_num()-1)
                osi[1].append([si])
                sio[1][si].append(len(osi[0])-1)
                pointer.append(len(osi[0])-1)

        # Adjustment helper function
        def adjust(ox, dist=None):
            # Process input
            if not isinstance(ox, list):
                ox = [ox]

            for o in ox:
                # Get positions
                posSi = mol.pos(sio[0][osi[1][o][0]])
                posO = mol.pos(osi[0][o])

                # Determine new position
                if dist is not None:
                    pos = posSi
                    pos[2] += dist
                else:
                    pos = posO
                    pos[1] = posSi[1]

                # Move atom
                mol.put(osi[0][o], pos)

        # Run through silicon molecules
        for i in range(len(sio[0])):
            # x-axis
            if drill == "x":
                if abs(data[2][sio[0][i]]-box[2]) < 10e-3:
                    for o in sio[1][i]:
                        if abs(data[2][osi[0][o]]-gap[2]) < 10e-3:
                            add_oxy(i, o, si_grid)
                            adjust(o, -si_grid)

            # y-axis
            elif drill == "y":
                # Silica at the front
                if abs(data[2][sio[0][i]]-box[2]+y_grid) < 10e-3:
                    for o in sio[1][i]:
                        if abs(data[2][osi[0][o]]-gap[2]) < 10e-3:
                            add_oxy(i, o, si_grid)
                            adjust(o, -si_grid)

                # Silica at the back
                if abs(data[2][sio[0][i]]-gap[2]) < 10e-3:
                    ox = []
                    theta = []
                    for o in sio[1][i]:
                        # Single oxygen
                        if abs(data[2][osi[0][o]]-box[2]+y_grid) < 10e-3:
                            add_oxy(i, o, -si_grid)
                            adjust(o, si_grid)
                            break

                        # Double oxygens
                        if abs(data[2][osi[0][o]]-box[2]) < 10e-3:
                            # Find second silicon and get positions
                            si = osi[1][o][0] if not osi[1][o][0] == i else osi[1][o][1]
                            oP = data[0][osi[0][o]]
                            sP = data[0][sio[0][si]]

                            # Add to list
                            ox.append(o)

                            # Left, boundary left and right oxygen (this order)
                            theta.append(angle-180) if sP > oP or abs(sP -
                                                                      oP) > box[2]/2 else theta.append(180-angle)

                    add_oxy(i, ox, -si_grid, theta)
                    adjust(ox)

            # z-axis
            elif drill == "z":
                if abs(data[2][sio[0][i]]-gap[2]) < 10e-3:
                    for o in sio[1][i]:
                        if abs(data[2][osi[0][o]]-box[2]) < 10e-3:
                            add_oxy(i, o, -si_grid)

    def drill(self, focal, diam):
        """Drill through a pore through the silicon-oxygen-grid. This function
        does not delete any atoms but rather unlinks all bonds of atoms
        within the pore diameter.

        Parameters
        ----------
        focal : list
            Focal point of the grid
        diam : float
            Pore diameter
        """
        # Initialize
        data = self._data
        dim = self._dim
        sio = self._sio
        osi = self._osi

        # Delete silicon atoms
        lists = [sio, osi]

        for i in range(len(lists)):
            for target in lists[i][0]:
                length = math.sqrt(sum([(data[x][target]-focal[x])**2 for x in range(dim-1)]))

                if length < diam/2:
                    self.unlink(target)

    def prepare(self):
        """Prepare binding site in the pore. These are oxygen atoms with only
        one bond to a grid silicon and another loose end directed to the pore.

        Inspired by earlier work, the pore is then prepared in the steps

        1. Remove all unsaturated silicon, meaning all silicon with less than four bonds
        2. Remove all silicon with three binding sites

        Thus only silicon are left with a maximum of two binding sites.
        """
        # Initialize
        sio = self._sio
        osi = self._osi

        # Remove unsaturated silicon
        for i in range(len(sio[0])):
            if len(sio[1][i]) < 4:
                self.unlink(sio[0][i])

        # Remove silicon with three binding sites (2 times for safety)
        for i in range(2):
            temp_o = []
            for i in range(len(osi[0])):
                if len(osi[1][i]) == 1:
                    temp_o.extend(osi[1][i])

            temp_s = list(set(temp_o))

            for s in temp_s:
                if temp_o.count(s) > 2:
                    self.unlink(sio[0][s])

    def site(self):
        """Find and return the matrix :math:`\\boldsymbol{B}` of binding sites

        .. math::

            \\boldsymbol{B}=\\begin{bmatrix}
                o_0&s_{0,0}&t_0&g_0\\\\
                o_1&s_{1,0}&t_1&g_1\\\\
                \\vdots\\\\
                o_b&s_{b,0}&t_b&g_b
            \\end{bmatrix}

        with number of binding sites :math:`b`. The type
        :math:`t_{k=0,\\dots,b}` is 0 for binding sites inside the pore and 1
        for sites on the outside. Geminal entries :math:`g_k` contain pointers
        to binding sites in the binding matrix :math:`\\boldsymbol{B}`.

        Returns
        -------
        site : list
            Binding site matrix
        """
        # Initialize
        osi = self._osi
        sio = self._sio
        data = self._data
        pointer = self._pointer
        geminal = self._geminal()
        box = self._mol.get_box_c()
        si_grid = self._si_grid
        site = []

        # Find binding sites - 0-O, 1-Si, 2-type, 3-geminal
        for i in range(len(osi[0])):
            if len(osi[1][i]) == 1:
                # Add indices
                entry = [osi[0][i], sio[0][osi[1][i][0]]]

                # Add type
                zCoord = data[2][entry[0]]
                if (abs(zCoord) < 10e-3 or abs(zCoord-box[2]) < 10e-3 or
                    abs(zCoord-si_grid) < 10e-3 or abs(zCoord-box[2]+si_grid) < 10e-3 or
                        abs(zCoord-0.107) < 10e-3 or abs(zCoord-box[2]+0.082) < 10e-3):
                    entry.append(1)
                else:
                    entry.append(0)

                # Add geminal
                entry.append(None)

                # Add entry to site
                site.append(entry)

        # Fill geminal
        temp_site = utils.column(site)
        for gem in geminal:
            indices = [i for i, x in enumerate(temp_site[1]) if x == gem]
            site[indices[0]][3] = indices[1]
            site[indices[1]][3] = indices[0]

        return site
