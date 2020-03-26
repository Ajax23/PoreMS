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

from porems.cube import Cube
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
        self._pore = self.get_pattern().generate(self._size, self._drill)
        self._pore.set_name("pore")
        self._t_tot["Block"] = utils.toc(t, "Block   ", is_time)

        # Cube and bonding
        t = utils.tic()
        self._cube = Cube(self._pore, cube_size, True) # Create cubes
        #self._bonding = Bonding(self._verlet)          # Create bond matrix
        self._t_tot["Matrix"] = utils.toc(t, "Matrix  ", is_time)


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


    ##################
    # Getter Methods #
    ##################
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

    def get_pore(self):
        """Return the finished pore.

        Returns
        -------
        pore : Molecule
            Finished pore object
        """
        return self._pore
