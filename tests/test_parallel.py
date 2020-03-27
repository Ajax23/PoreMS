import os
import sys

import unittest

from porems.atom import Atom
from porems.molecule import Molecule
from porems.essentials import *
from porems.store import Store
from porems.pattern import *
from porems.cube import Cube
from porems.matrix import Matrix
from porems.pore import Pore


class UserModelCase(unittest.TestCase):
    ########
    # Cube #
    ########
    def test_cube_parallel(self):
        cube = Cube(BetaCristobalit().generate([2, 2, 2], "z"), 0.4, True)

        self.assertEqual(len(cube.find_parallel(None, ["Si", "O"], 0.155, 0.005)), 192)
        self.assertEqual(len(cube.find_parallel(None, ["O", "Si"], 0.155, 0.005)), 384)


if __name__ == '__main__':
    unittest.main(verbosity=2)
