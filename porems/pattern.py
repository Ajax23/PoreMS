################################################################################
# Pattern Pack                                                                 #
#                                                                              #
"""This file contains minimal pattern structures and their properties for
creating pore structures."""
################################################################################


import math
import copy

from porems.molecule import Molecule


class Pattern():
    """This class is a container for individual pattern classes."""
    def __init__(self):
        self._repeat = [0, 0, 0]
        self._gap = [0, 0, 0]

    def pattern(self):
        """Construct minimal block structure.

        Returns
        -------
        block : Molecule
            Minimal block structure
        """
        return Molecule()

    def get_repeat(self):
        """Return the repeat distances in all dimensions.

        Returns
        -------
        repeat : list
            List of repetition distances
        """
        return self._repeat

    def get_gap(self):
        """Return the gap between the pore block and the box edge in all
        dimensions.

        Returns
        -------
        gap : list
            List of edge distances
        """
        return self._gap


class BetaCristobalit(Pattern):
    """This class defines the minimal structure of a :math:`\\beta`-crystobalit
    molecule.
    """
    def __init__(self):
        # Call super class
        super(BetaCristobalit, self).__init__()

        # Set bond length si-o and bond angle si-o-si
        self._b = 0.155
        self._a = 2*math.atan(math.sqrt(2))*180/math.pi

        # Set repetition distances and gap towards the box edge
        self._repeat = [0.506, 0.877, 1.240]
        self._gap = [0.126, 0.073, 0.155]

    def _hexagonal(self):
        """Define hexagonal molecule.

        Returns
        -------
        hex : Molecule
            Hexagonal molecule object
        """
        hex = Molecule()

        hex.add("Si", [0, 0, 0])
        hex.add("O",   0, r=self._b, theta=self._a, phi=60)
        hex.add("Si",  1, r=self._b, bond=[0, 1])
        hex.add("O",   2, r=self._b, theta=180,  phi=0)
        hex.add("Si",  3, r=self._b, bond=[2, 3])
        hex.add("O",   4, r=self._b, theta=self._a, phi=300)
        hex.add("Si",  5, r=self._b, bond=[4, 5])
        hex.add("O",   6, r=-self._b, theta=self._a, phi=60)
        hex.add("Si",  7, r=self._b, bond=[6, 7])
        hex.add("O",   8, r=-self._b, theta=180,  phi=0)
        hex.add("Si",  9, r=self._b, bond=[8, 9])
        hex.add("O",  10, r=-self._b, theta=self._a, phi=300)

        hex.rotate("y", 90)
        hex.rotate("z", 90)
        hex.rotate("y", 180)
        hex.rotate("x", hex.angle(hex.bond(8, 6)[0], [0, 1, 0]))
        hex.rotate("x", -90)

        return hex

    def pattern(self):
        """Construct minimal block structure.

        Returns
        -------
        block : Molecule
            Minimal block structure
        """
        # Initialize
        mols = [self._hexagonal() for x in range(3)]

        # Combine three hexagonal molecules
        mols[0].add("O",  2, r=self._b)
        mols[0].add("O",  6, r=self._b)
        mols[0].add("O", 10, r=self._b)

        mols[1].move(0, mols[0].pos(12))
        mols[1].translate([0, 0, self._b])
        mols[1].delete([1, 2, 3, 4, 5, 6, 7])

        mols[2].move(0, mols[0].pos(14))
        mols[2].translate([0, 0, self._b])
        mols[2].delete([2, 3, 4, 5, 6, 7, 8, 9, 10, 11])

        # Build block
        block = []
        for i in range(11):
            block.append(Molecule(inp=mols))

        block[1].rotate("x", 180)
        block[1].move(18, block[0].pos(18))
        block[1].translate([0, 0, self._b*2])
        block[1].add("O", block[0].pos(18), r=self._b)

        block[2].rotate("x", 180)
        block[2].move(0, block[0].pos(18))

        block[3].move(4, block[0].pos(15))

        block[4].rotate("x", 180)
        block[4].move(16, block[0].pos(18))
        block[4].delete([21, 19, 15, 20, 14, 12, 10, 2, 11, 1, 0])

        block[5].move(6, block[1].pos(16))
        block[5].delete([21, 19, 15, 20, 14, 12, 10, 2, 11, 1, 0])

        # Delete overlaping molecules for repetition
        block = Molecule(inp=block)
        block.overlap()

        # Determine number of needed blocks
        block.zero()

        # Delete overlapping atoms
        block.delete([4, 3, 2, 12, 15, 45, 46, 55, 57, 63, 26, 25, 24, 34, 37, 61])  # x-Axis
        block.delete([50, 48, 36, 40, 41])  # y-Axis
        block.delete([45, 20, 19, 21, 22, 23, 24, 25, 17, 18])  # z-Axis

        # Move to zero
        block.zero()

        return block
