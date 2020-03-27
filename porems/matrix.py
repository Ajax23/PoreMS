################################################################################
# Matric Class                                                                 #
#                                                                              #
"""Matrix class for preserving bond information."""
################################################################################


import math

import porems.utils as utils


class Matrix:
    """The aim of this class is preserving all information of the grid bonds.

    Although the search can be parallelized, still multiple iterations are
    needed to cover the surface preparations. Additionally, due to machine
    inaccuracy there is the risk of bonds not being detected as such, leading
    to artefacts. Also it is not possible to ensure that all bonds were  found
    when deleting atoms because all systems are shaped differently.  Therefore,
    another optimization, or rather supporting algorithm, was implemented to
    bypass these issues.

    The idea was reducing the number of iterations to a single search by
    creating a connectivity matrix of all oxygen and silicon atoms. Therefore
    two matrices :math:`\\boldsymbol{O}` for oxygen and :math:`\\boldsymbol{S}`
    for silicon were defined

    .. math::

        \\begin{array}{cc}
            \\boldsymbol{O}=
            \\begin{bmatrix}
                o_0&\\begin{pmatrix}p_{s,0,0}&p_{s,0,1}\\end{pmatrix}\\\\
                o_1&\\begin{pmatrix}p_{s,1,0}&p_{s,1,1}\\end{pmatrix}\\\\
                \\vdots&\\vdots\\\\
                o_k&\\begin{pmatrix}p_{s,k,0}&p_{s,k,1}\\end{pmatrix}
            \\end{bmatrix}
            ,&
            \\boldsymbol{S}=
            \\begin{bmatrix}
                s_0&\\begin{pmatrix}p_{o,0,0}&p_{o,0,1}&p_{o,0,2}&p_{o,0,3}\\end{pmatrix}\\\\
                s_1&\\begin{pmatrix}p_{o,1,0}&p_{o,1,1}&p_{o,1,2}&p_{o,1,3}\\end{pmatrix}\\\\
                \\vdots&\\vdots\\\\
                s_l&\\begin{pmatrix}p_{o,l,0}&p_{o,l,1}&p_{o,l,2}&p_{o,l,3}\\end{pmatrix}
            \\end{bmatrix}
        \\end{array}

    with atom ids of oxygen :math:`o` and silicon :math:`s` in the data matrix
    of  the pore, list id pointer :math:`p_s` of the silicon entry in matrix
    :math:`\\boldsymbol{S}` and pointer :math:`p_o` of the oxygen entry in
    matrix :math:`\\boldsymbol{O}`. Thus each entry of the matrix presents the
    atom and all  its binding partners. These matrices are filled after creating
    the cristobalite  block. This way it is possible to check whether all bonds
    were found, since all  entries need to have the same size when considering
    periodic boundary conditions.

    Using this implementation, it is no longer required to physically delete
    atoms when carving out the pore, it is enough to remove binding partners
    from the matrices. Thus, the surface preparation conditions only need to
    consider the number of bonds remaining in each entry and thereby determine
    whether an atom needs to be removed or not, resulting into an effort scaling
    linear with the number of atoms

    .. math::

      \\mathcal{O}(n).

    Parameters
    ----------
    mol : Molecule
        Molecule grid object
    orient : string
        Orientation of the molecule block
    bonds : list
        List of all pairwise bonds in the grid with only two different atom
        types
    """
    def __init__(self, mol, orient, bonds):
        # Initialize
        self._dim = 3
        self._mol = mol

        # Create bond dictionary
        self._matrix = {}
        for bond in bonds:
            # Fill bond as given
            self._matrix[bond[0]] = bond[1]
            # Fill reverse bonds
            for atom_b in bond[1]:
                if not atom_b in self._matrix:
                    self._matrix[atom_b] = []
                self._matrix[atom_b].append(bond[0])


    ###################
    # Private Methods #
    ###################
    def _unbind(self, atom_a, atom_b):
        """Remove the bond from the connectivity matrix.

        Parameters
        ----------
        atom_a : integer
            Atom A
        atom_b : integer
            Atom B
        """
        self._matrix[atom_a].remove(atom_b)
        self._matrix[atom_b].remove(atom_a)

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
        y_grid = 0.073

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
                            theta.append(angle-180) if sP > oP or abs(sP - oP) > box[2]/2 else theta.append(180-angle)

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
