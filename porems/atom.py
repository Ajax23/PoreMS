################################################################################
# Atom Class                                                                   #
#                                                                              #
"""All necessary function for creating and editing atoms."""
################################################################################


import pandas as pd


class Atom:
    """This class defines an atom object containing the position, the atom type
    and optionally a specific name.

    Parameters
    ----------
    pos : list
        Atom position
    atom_type : string
        Atom type as in the periodic table of elements
    name : string, optional
        Atom name
    residue : integer, optional
        Residue number
    """
    def __init__(self, pos, atom_type, name="", residue=0):
        # Initialize
        self._pos = pos
        self._atom_type = atom_type
        self._name = name
        self._residue = residue


    ##################
    # Representation #
    ##################
    def __repr__(self):
        """Create a pandas table of the atom data.

        Returns
        -------
        repr : DataFrame
            Pandas data frame of the molecule object
        """
        # Set colums names
        columns = ["Residue", "Name", "Type", "x", "y", "z"]

        # Get data
        data =[[self._residue, self._name, self._atom_type,
                self._pos[0], self._pos[1], self._pos[2]]]

        # Create data frame
        return pd.DataFrame(data, columns=columns).to_string()


    ##################
    # Setter Methods #
    ##################
    def set_pos(self, pos):
        """Set the atom position.

        Parameters
        ----------
        pos : list
            Atom position
        """
        self._pos = pos

    def set_atom_type(self, atom_type):
        """Set the atom type.

        Parameters
        ----------
        atom_type : string
            Atom type as in the periodic table of elements
        """
        self._atom_type = atom_type

    def set_name(self, name):
        """Set the atom name.

        Parameters
        ----------
        name : string
            Atom name
        """
        self._name = name

    def set_residue(self, residue):
        """Set the residue index.

        Parameters
        ----------
        residue : integer
            Residue index
        """
        self._residue = residue


    ##################
    # Getter Methods #
    ##################
    def get_pos(self):
        """Return the atom position.

        Returns
        -------
        pos : list
            Atom position
        """
        return self._pos

    def get_atom_type(self):
        """Return the atom type.

        Returns
        -------
        atom_type : string
            Atom type
        """
        return self._atom_type

    def get_name(self):
        """Return the atom name.

        Returns
        -------
        name : string
            Atom name
        """
        return self._name

    def get_residue(self):
        """Return the residue index.

        Returns
        -------
        residue : integer
            Residue index
        """
        return self._residue
