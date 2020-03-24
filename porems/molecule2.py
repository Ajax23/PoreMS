################################################################################
# Molecule Class                                                               #
#                                                                              #
"""All necessary function for creating and editing molecules."""
################################################################################


import math
import pandas as pd

import porems.utils as utils
import porems.geometry as geometry

from porems.atom import Atom


class Molecule:
    """

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

    """
    def __init__(self, name="molecule", short="MOL", inp=None):
        # Initialize
        self._dim = 3

        self._name = name
        self._short = short

        self._box = None

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


    ################################
    # Private Methods - Management #
    ################################
    def _read(self, file_path, file_type):
        """Read a molecule from a file. Currently only **GRO**, **PDB** and
        **MOL2** files are supported.

        Parameters
        ----------
        file_path : string
            Link to requested file
        file_type : string
            Int for types list id or str for file extension name

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

                # Gro file
                if file_type == "GRO":
                    if line_idx > 0 and len(line_val) > 3:
                        pos = [float(line_val[i]) for i in range(3, 5+1)]
                        name = line_val[1]
                        atom_type = ''.join([i for i in line_val[1] if not i.isdigit()])

                        atom_list.append(Atom(pos, atom_type, name))

                # Pdb file
                elif file_type == "PDB":
                    if line_val[0] in ["ATOM", "HETATM"]:
                        pos = [float(line_val[i])/10 for i in range(6, 8+1)]
                        name = line_val[11]
                        atom_type = line_val[11]

                        atom_list.append(Atom(pos, atom_type, name))

                # Mol2 file
                elif file_type == "MOL2":
                    if len(line_val) > 8:
                        pos = [float(line_val[i])/10 for i in range(2, 4+1)]
                        name = line_val[1]
                        atom_type = ''.join([i for i in line_val[1] if not i.isdigit()])

                        atom_list.append(Atom(pos, atom_type, name))

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
            List of atoms to be duplicated

        Returns
        -------
        mol : Molecule
            Molecule object
        """
        return Molecule(inp=[self.get_atom_list[x] for x in atoms])

    def _append(self, mol):
        """Append a given molecule to the current object.

        Parameters
        ----------
        mol : Molecule
            Molecule object
        """
        self._atom_list += mol.get_atom_list

    def _column_pos(self):
        """Create column list of atom positions

        Returns
        -------
        column : list
            Columns of all atom positions in all dimensions
        """
        return utils.column([atom.get_pos() for atom in self._atom_list])


    ##############################
    # Private Methods - Geometry #
    ##############################
    def _vector(pos_a, pos_b):
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
        elif not isinstance(pos_a, list) and isinstance(pos_b, list):
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

            \\boldsymbol{b}=\\begin{pmatrix}\\max(\\boldsymbol{d}_1)&max(\\boldsymbol{d}_1)&\\dots&max(\\boldsymbol{d}_n)\\end{pmatrix}^T

        where :math:`\\boldsymbol{d}_i` is the dimension-vector of the data
        matrix.

        Returns
        -------
        box : list
            Box length of the current molecule
        """
        data = self._column_pos()
        return [max(data[i]) if max(data[i]) > 0 else 0.001 for i in range(self._dim)]


    ##############################
    # Public Methods - Transform #
    ##############################



    #########################
    # Public Methods - Edit #
    #########################



    ##############################
    # Public Methods - Calculate #
    ##############################
    def centroid(self):
        """Calculate the geometrical center of mass

        .. math::

            \\text{com}=\\begin{pmatrix}c_1&c_2&\\dots&c_n\\end{pmatrix}^T

        with

        .. math::

            c_i=\\frac{1}{m}\\sum_j^m d_{ij}

        Returns
        -------
        Centroid : list
            Geometrical center of mass
        """
        # Calculate the centroid
        data = self._column_pos()
        return [sum(data[i])/len(data[i]) for i in range(self._dim)]



    #########################
    # Public Methods - Edit #
    #########################



    ##################
    # Public Methods #
    ##################
    def __repr__(self):
        """Create a pandas table of the molecule data."""
        # Create data table from atom list
        atom_data = []
        for atom in self._atom_list:
            atom_data.append([atom.get_name(), atom.get_atom_type()])
            atom_data[-1].extend([atom.get_pos()[i] for i in range(self._dim)])
        atom_data = utils.column(atom_data)

        # Create column names
        column_names = ["Name", "Type"]
        column_names.extend([["x", "y", "z"][dim] for dim in range(self._dim)])

        # Create dictionary
        data = {column_names[i]: atom_data[i] for i in range(self._dim+2)}

        # Create dataframe
        return pd.DataFrame(data)

    def __str__(self):
        """Convert a pandas table of the molecule data to string format."""
        return self.__repr__().to_string()


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


    ##################
    # Getter Methods #
    ##################
    def get_atom_list(self):
        """Return the atoms list.

        Returns
        -------
        atom_list : list
            Atom list
        """
        return self._atom_list

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
        """Return the boxsize of the molecule.

        Returns
        -------
        box : list
            Box size in all dimensions
        """
        return self._box_size() if self._box is None else self._box
