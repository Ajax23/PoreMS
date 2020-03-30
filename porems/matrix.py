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
    to artefacts. Also it is not possible to ensure that all bonds were found
    when deleting atoms, because all systems are shaped differently. Therefore,
    another optimization, or rather supporting algorithm, was implemented to
    bypass these issues.

    The idea was reducing the number of iterations to a single search by
    creating a connectivity matrix of all grid atoms. The result is a dictionary
    that has atoms :math:`1\\dots n` as keys and their corresponding value is a
    list of bonded atoms :math:`1\\dots m`

    .. math::

        \\boldsymbol{C}=
        \\begin{bmatrix}
            \\text{atom}_1&\\begin{pmatrix}\\text{atom}_{1,1}&\\text{atom}_{1,2}&\\dots&\\text{atom}_{1,m_1}\\end{pmatrix}\\\\
            \\text{atom}_2&\\begin{pmatrix}\\text{atom}_{2,1}&\\text{atom}_{2,2}&\\dots&\\text{atom}_{2,m_2}\\end{pmatrix}\\\\
            \\vdots&\\vdots\\\\
            \\text{atom}_n&\\begin{pmatrix}\\text{atom}_{n,1}&\\text{atom}_{n,2}&\\dots&\\text{atom}_{n,m_n}\\end{pmatrix}\\\\
        \\end{bmatrix}

    Using this implementation, it is no longer required to physically delete
    atoms when carving out a structure, it is enough to remove binding partners
    from the matrix. For example conditions for the surface preparation only
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
            List atoms or one atom id
        """
        # Porocess input
        atoms = [atoms] if isinstance(atoms, int) else atoms

        # Split all bonds
        for atom_a in atoms:
            atoms_b = self._matrix[atom_a][:]
            for atom_b in atoms_b:
                self. split(atom_a, atom_b)

    def bound(self, num_bonds):
        """Return a list of atoms with the specified number of bonds.

        Parameters
        ----------
        num_bonds : int
            Number of bonds to search for
        """
        return [atom for atom in self._matrix if len(self._matrix[atom])==num_bonds]


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