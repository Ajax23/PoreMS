################################################################################
# Pore Class                                                                   #
#                                                                              #
"""Extension of the molecule class for pores."""
################################################################################


import math
import copy
import random
import multiprocessing as mp

import porems.utils as utils

from porems.molecule import Molecule
from porems.pattern import BetaCristobalit

from porems.dice import Dice
from porems.matrix import Matrix
from porems.shape import Cylinder


class Pore():
    """Pore class

    Parameters
    ----------
    size : list
    orient : string

    Examples
    --------

    """
    def __init__(self, size, orient):
        # Initialize
        self._dim = 3
        self._size = size
        self._orient = orient

        self._pattern = None
        self._shape = None
        self._t_tot = {}


    ############
    # Building #
    ############



    #################
    # Binding sites #
    #################



    #######################
    # Molecule Attachment #
    #######################



    ###############
    # Final Edits #
    ###############



    ###########
    # Analyze #
    ###########



    ##############
    # Generation #
    ##############
    def generate(self, cube_size=0.4, is_time=True):
        """Run pore generation.

        Parameters
        ----------
        cube_size : float, optional
            Cube size
        is_time : bool
            True to print time usage
        """
        # Create block structure
        t = utils.tic()
        self._pore = self.get_pattern().generate(self._size, self._orient)
        self._pore.set_name("pore")
        self._t_tot["Block"] = utils.toc(t, "Block   ", is_time)

        # Cube and bonding
        # t = utils.tic()
        # dice = Dice(self._pore, cube_size, True)
        # bonds = dice.find_parallel(None, ["Si", "O"], 0.155, 10e-2)
        # matrix = Matrix(bonds)
        # shape = Cylinder(self._pore.centroid()[:2]+[0], [1, 1, 1], {"length": 6, "diameter": 4})
        #
        # shape.plot()
        #
        # for atom_id, atom in enumerate(self._pore.get_atom_list()):
        #     if shape.is_in(atom.get_pos()):
        #         matrix.strip(atom_id)
        #
        # self._pore.delete(matrix.bound(0))

        # # Consistency check
        # if not len(bond[1])==atom_type[0][1]:
        #     print("Inconsistency in number of bonds for "+str(bond)+"...")
        #     return


        #self._bonding = Bonding(self._verlet)          # Create bond matrix
        # self._t_tot["Matrix"] = utils.toc(t, "Matrix  ", is_time)


    ##################
    # Setter Methods #
    ##################
    def set_pattern(self, pattern=None):
        """Set the block pattern.

        Parameters
        ----------
        pattern : Pattern, None, optional
            Pattern object
        """
        self._pattern = BetaCristobalit() if pattern is None else pattern

    def set_shape(self, shape=None):
        """Set the pore shape.

        Parameters
        ----------
        shape : Shape, None, optional
            Shape object
        """
        self._shape = shape


    ##################
    # Getter Methods #
    ##################
    def get_pore(self):
        """Return the finished pore.

        Returns
        -------
        pore : Molecule
            Finished pore object
        """
        return self._pore

    def get_pattern(self):
        """Return the pattern object.

        Returns
        -------
        pattern : Pattern
            Pattern object
        """
        if self._pattern is None:
            self.set_pattern()
        return self._pattern

    def get_shape(self):
        """Return the shape object.

        Returns
        -------
        shape : Shape
            shape object
        """
        return self._shape
