################################################################################
# Pattern Pack                                                                 #
#                                                                              #
"""This file contains minimal pattern structures and their properties for
creating pore structures."""
################################################################################


import math
import copy

import porems.geometry as geometry

from porems.molecule import Molecule


class Pattern():
    """This class is a container for individual pattern classes."""
    def __init__(self):
        self._dim = 3
        self._repeat = [0, 0, 0]
        self._gap = [0, 0, 0]
        self._size = [0, 0, 0]


    ##############
    # Generation #
    ##############
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

            for i in range(self._num[dim]):
                temp = copy.deepcopy(block)
                vec = [(i+1)*self._repeat[dim] if j == dim else 0 for j in range(self._dim)]
                temp.translate(vec)
                p.append(temp)
            self._block(dim+1, Molecule(inp=p))
        else:
            self._structure = block

    def _orientation(self, orient):
        """Rotate pore orientation, so that the specified axis becomes the
        z-axis.

        Parameters
        ----------
        orient : string
            Axis that will be oriented towards the z-axis
        """
        if orient == "x":
            self._structure.rotate("y", 90)
            self._gap[0], self._gap[2] = self._gap[2], self._gap[0]
            self._size[0], self._size[2] = self._size[2], self._size[0]
        elif orient == "y":
            self._structure.rotate("x", 90)
            self._gap[1], self._gap[2] = self._gap[2], self._gap[1]
            self._size[1], self._size[2] = self._size[2], self._size[1]

        self._orient = orient

    def generate(self, size, orient):
        """Generate full block structure.

        Parameters
        ----------
        size : list
            Desired block dimensions
        orient : string
            Orientation of the block that should be drilled

        Returns
        -------
        structure : Molecule
            Full block structure
        """
        # Calculate repetition and size
        self._num = [round(size[i]/self._repeat[i]) for i in range(self._dim)]
        self._size = [self._repeat[i]*self._num[i] for i in range(self._dim)]

        # Generate block
        self._block(0, self.pattern())

        # Rotate block towards desired orientation
        self._orientation(orient)

        # Translate gap
        self._structure.zero()
        self._structure.translate([x/2 for x in self._gap])
        self._structure.set_box(self._size)

        return self._structure


    ##################
    # Getter Methods #
    ##################
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

    def get_size(self):
        """Return the size of the generated block structure.

        Returns
        -------
        size : list
            List of molecule dimentsions
        """
        return self._size

    def get_block(self):
        """Return the generated molecule object.

        Returns
        -------
        block : Molecule
            Generated block molecule
        """
        return self._structure

    def get_orient(self):
        """Return the orientation of the block.

        Returns
        -------
        orient : string
            Block orientation
        """
        return self._orient


class BetaCristobalit(Pattern):
    """This class defines the minimal structure of a :math:`\\beta`-cristobalite
    molecule."""
    def __init__(self):
        # Call super class
        super(BetaCristobalit, self).__init__()

        # Set bond length si-o and bond angle si-o-si
        self._b = 0.155
        self._a = 2*math.atan(math.sqrt(2))*180/math.pi

        # Set repetition distances and gap towards the box edge
        self._repeat = [0.506, 0.877, 1.240]
        self._gap = [0.126, 0.073, 0.155]


    ############
    # Building #
    ############
    def _hexagonal(self):
        """Define hexagonal molecule.

        Returns
        -------
        hex : Molecule
            Hexagonal molecule object
        """
        hex = Molecule()

        hex.add("Si", [0, 0, 0])
        hex.add("O",   0, r= self._b, theta=self._a, phi=60)
        hex.add("Si",  1, r= self._b, bond=[0, 1])
        hex.add("O",   2, r= self._b, theta=180,  phi=0)
        hex.add("Si",  3, r= self._b, bond=[2, 3])
        hex.add("O",   4, r= self._b, theta=self._a, phi=300)
        hex.add("Si",  5, r= self._b, bond=[4, 5])
        hex.add("O",   6, r=-self._b, theta=self._a, phi=60)
        hex.add("Si",  7, r= self._b, bond=[6, 7])
        hex.add("O",   8, r=-self._b, theta=180,  phi=0)
        hex.add("Si",  9, r= self._b, bond=[8, 9])
        hex.add("O",  10, r=-self._b, theta=self._a, phi=300)

        hex.rotate("y", 90)
        hex.rotate("z", 90)
        hex.rotate("y", 180)
        hex.rotate("x", geometry.angle(hex.bond(8, 6), [0, 1, 0]))
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
            block.append(Molecule(inp=copy.deepcopy(mols)))

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

        # Delete overlapping molecules
        block = Molecule(inp=block)
        overlap = block.overlap()
        block.delete(sum([overlap[x] for x in overlap], []))

        # Check overlapping atoms due to repetition
        for dim in range(3):
            for pm in range(2):
                # Define translate vector
                translate = [0, 0, 0]
                translate[dim] = self._repeat[dim]
                translate[dim] *= -1 if pm == 1 else 1

                # Remove atoms
                block = Molecule(inp=[block])
                mol_repeat = copy.deepcopy(block)
                mol_repeat.translate(translate)
                block.delete([x for x in  Molecule(inp=[block, mol_repeat]).overlap() if x < block.get_num()])

        # Move to zero
        block.zero()

        return block
