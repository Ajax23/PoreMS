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

from porems.pattern import *

from porems.molecule import Molecule
from porems.verlet import Verlet
from porems.bonding import Bonding


class Pore():
    """

    Parameters
    ----------

    Examples
    --------

    """
    def __init__(self, size, drill):
        # Initialize
        self._dim = 3
        self._size = size if isinstance(size, list) else [size for s in range(self._dim)]
        self._drill = drill

        self._pattern = None
        self._t_tot = {}


    ############
    # Building #
    ############
    def _block(self, dim, block):
        """Recursively duplicate and translate a given molecule block in all
        given dimensions. The duplication stops if the next added block would
        create the pore longer than specified in the constructor.

        Parameters
        ----------
        dim : integer
            Repeat dimensions
        block : Molecule
            Molecule unit to be duplicated
        """
        if dim < self._dim:
            p = []

            for i in range(self._num_rep[dim]):
                temp = copy.deepcopy(block)
                vec = [(i+1)*self._repeat[dim] if j == dim else 0 for j in range(self._dim)]

                temp.translate(vec)

                p.append(temp)

            self._block(dim+1, Molecule(inp=p))
        else:
            self._pore = block

    def _orientation(self, pore):
        """Rotate pore orientation, so that the specified drill axis becomes
        the z-axis.

        Parameters
        ----------
        pore : Molecule
            Pore object to rotate
        """
        # Initialize
        drill = self._drill
        gap = self._gap
        size = self._size

        # Rotate pore
        if drill == "x":
            pore.rotate("y", 90)
        elif drill == "y":
            pore.rotate("x", 90)

        # Set zero
        pore.zero()

        # Update gap and size lists
        if drill == "x":
            self._gap = [gap[2], gap[1], gap[0]]
            self._size = [size[2], size[1], size[0]]
        elif drill == "y":
            self._gap = [gap[0], gap[2], gap[1]]
            self._size = [size[0], size[2], size[1]]


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
    # Analyse #
    ###########



    ################
    # Finalization #
    ################
    def generate(self):
        # Set pattern
        if self._pattern is None:
            self.set_pattern(BetaCristobalit())

        # Create block structure
        t = utils.tic()
        self._block(0, self._pattern)                  # Build block
        self._orientation(self._pore)                  # Rotate drill axis
        self._pore.translate(self._gap)                # Translate gap
        self._pore.set_name("pore")                    # Set pore name
        self._t_tot["Build"] = utils.toc(t, "Build   ", False)

    ##################
    # Setter Methods #
    ##################
    def set_pattern(self, pattern):
        """Set the block pattern.

        Parameters
        ----------
        pattern : Pattern
            Pattern object
        """
        # Get information
        self._repeat = pattern.get_repeat()
        self._gap = pattern.get_gap()
        self._pattern = pattern.pattern()

        # Calculate repetition and size
        self._num_rep = [round(self._size[i]/self._repeat[i]) for i in range(self._dim)]
        self._size = [self._repeat[i]*self._num_rep[i]+self._gap[i] for i in range(self._dim)]


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
