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


    ############################
    # Public Methods - Bonding #
    ############################
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
        self._dice = Dice(self._pore, cube_size, True) # Create cubes

        # # Consistency check
        # if not len(bond[1])==atom_type[0][1]:
        #     print("Inconsistency in number of bonds for "+str(bond)+"...")
        #     return


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
