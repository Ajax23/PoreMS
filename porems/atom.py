################################################################################
# Atom Class                                                                   #
#                                                                              #
"""All necessary function for creating and editing atoms."""
################################################################################


import pandas as pd


class Atom:
    """This class defines an atom object, containing position and atom_type.

    Parameters
    ----------
    pos : list
        Atom position
    atom_type : string
        Atom type as in PEM
    name : string, None, optional
        Atom name

    Examples
    --------
    """
    def __init__(self, pos, atom_type, name=None):
        # Initialize
        self._pos = pos
        self._atom_type = atom_type
        self._name = name


    ##################
    # Representation #
    ##################
    def __repr__(self):
        """Create a pandas table of the molecule data."""
        # Create data table from atom list
        atom_data = [[self._name], [self._atom_type]]
        atom_data.extend([[self._pos[i]] for i in range(len(self._pos))])

        # Create column names
        column_names = ["Name", "Type"]
        column_names.extend([["x", "y", "z"][dim] for dim in range(len(self._pos))])

        # Create dictionary
        data = {column_names[i]: atom_data[i] for i in range(len(self._pos)+2)}

        # Create dataframe
        return pd.DataFrame(data)

    def __str__(self):
        """Convert a pandas table of the molecule data to string format."""
        return self.__repr__().to_string()


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
            Atom type as in PEM
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
        """Return the atom position.

        Returns
        -------
        name : string
            Atom name
        """
        return self._name
