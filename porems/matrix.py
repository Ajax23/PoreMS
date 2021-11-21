################################################################################
# Matric Class                                                                 #
#                                                                              #
"""Matrix class for preserving bond information."""
################################################################################


class Matrix:
    """The aim of this class is preserving all information of the grid bonds.

    Although the search can be parallelized, still multiple iterations are
    needed to cover the surface preparations. Additionally, due to machine
    inaccuracy there is the risk of bonds not being detected as such, leading
    to artefacts. Also, it is not possible to ensure that all bonds were found
    when deleting atoms, because all systems are shaped differently. Therefore,
    another optimization, or rather supporting algorithm, was implemented to
    bypass these issues.

    The idea was reducing the number of iterations to a single search by
    creating a connectivity matrix of all grid atoms. The result is a dictionary
    that has atoms :math:`1\\dots n` as keys and their corresponding value is a
    list of bonded atoms :math:`1\\dots m`

    .. math::

        \\boldsymbol{C}=
        \\begin{Bmatrix}
            a_1:&\\begin{bmatrix}a_{1,1}&a_{1,2}&\\dots&a_{1,m_1}\\end{bmatrix}\\\\
            a_2:&\\begin{bmatrix}a_{2,1}&a_{2,2}&\\dots&a_{2,m_2}\\end{bmatrix}\\\\
            \\vdots&\\vdots\\\\
            a_n:&\\begin{bmatrix}a_{n,1}&a_{n,2}&\\dots&a_{n,m_n}\\end{bmatrix}\\\\
        \\end{Bmatrix}

    Using this implementation, it is no longer required to physically delete
    atoms when carving out a structure, it is enough to remove binding partners
    from the matrix. For example, conditions for the surface preparation only
    need to consider the number of bonds remaining in each entry and thereby
    determine whether an atom needs to be removed or not, resulting into a
    negligible computational effort scaling linear with the number of atoms

    .. math::

      \\mathcal{O}(n).

    Parameters
    ----------
    bonds : list
        List of all pairwise bonds in the grid with only two different atom
        types
    """
    def __init__(self, bonds):
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


    ###########
    # Editing #
    ###########
    def split(self, atom_a, atom_b):
        """Remove the bond between two atoms from the connectivity matrix.

        Parameters
        ----------
        atom_a : integer
            Atom A
        atom_b : integer
            Atom B
        """
        self._matrix[atom_a].remove(atom_b)
        self._matrix[atom_b].remove(atom_a)

    def strip(self, atoms):
        """Remove all bonds of a specified atom from the connection matrix.

        Parameters
        ----------
        atoms : list, integer
            List of atom ids or one atom id
        """
        # Porocess input
        atoms = [atoms] if isinstance(atoms, int) else atoms

        # Split all bonds
        for atom_a in atoms:
            atoms_b = self._matrix[atom_a][:]
            for atom_b in atoms_b:
                self. split(atom_a, atom_b)

    def add(self, atom_a, atom_b):
        """Add atom between given atom ids.

        Parameters
        ----------
        atom_a : int
            First atom id
        atom_b : int
            Second atom id
        """
        # Add entry fort first atom
        if atom_a in self._matrix.keys():
            self._matrix[atom_a].append(atom_b)
        else:
            self._matrix[atom_a] = [atom_b]

        # Add entry for second atom
        if atom_b in self._matrix.keys():
            self._matrix[atom_b].append(atom_a)
        else:
            self._matrix[atom_b] = [atom_a]

    def bound(self, num_bonds, logic="eq"):
        """Return a list of atoms with the specified number of bonds. Possible
        ``logic`` arguments are

        * **eq** - Equals
        * **lt** - Less than
        * **gt** - Greater than

        Parameters
        ----------
        num_bonds : int
            Number of bonds to search for
        logic : string, optional
            Logic statement

        Returns
        -------
        atoms : list
            List of atom ids with the number of specified bonds
        """
        if logic=="eq":
            return [atom for atom in self._matrix if len(self._matrix[atom])==num_bonds]
        elif logic=="lt":
            return [atom for atom in self._matrix if len(self._matrix[atom])<num_bonds]
        elif logic=="gt":
            return [atom for atom in self._matrix if len(self._matrix[atom])>num_bonds]
        else:
            print("Matrix: Wrong logic statement...")
            return


    ##################
    # Getter Methods #
    ##################
    def get_matrix(self):
        """Return the connectivity matrix.

        Returns
        -------
        matrix : dictionary
            Connectivity matrix
        """
        return self._matrix
