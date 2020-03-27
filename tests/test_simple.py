import os
import sys

import unittest

# Install package
if sys.platform == "win32":
    os.system("pip install ../.")
else:
    os.system("pip install ../. &> /dev/null")
print("Finished inistalling package...")

# Import package
from porems.atom import Atom
from porems.molecule import Molecule
from porems.essentials import *
from porems.store import Store
from porems.pattern import *
from porems.cube import Cube
from porems.pore import Pore


class UserModelCase(unittest.TestCase):
    ########
    # Atom #
    ########
    def test_atom(self):
        atom = Atom([0.0, 0.1, 0.2], "H", "HO1")

        self.assertEqual(atom.__str__(), "  Name Type    x    y    z\n0  HO1    H  0.0  0.1  0.2")


    ############
    # Molecule #
    ############
    def test_molecule_loading(self):
        mol_gro = Molecule(inp="data/benzene.gro")
        mol_pdb = Molecule(inp="data/benzene.pdb")
        mol_mol2 = Molecule(inp="data/benzene.mol2")

        mol_atom = Molecule(inp=mol_mol2.get_atom_list())
        mol_concat = Molecule(inp=[mol_gro, mol_pdb])

        mol_append = Molecule(inp="data/benzene.gro")
        mol_append.append(mol_gro)

        pos_gro = [[round(x, 4) for x in col] for col in mol_gro.column_pos()]
        pos_pdb = [[round(x, 4) for x in col] for col in mol_pdb.column_pos()]
        pos_mol2 = [[round(x, 4) for x in col] for col in mol_mol2.column_pos()]
        pos_atom = [[round(x, 4) for x in col] for col in mol_atom.column_pos()]
        pos_concat = [[round(x, 4) for x in col] for col in mol_concat.column_pos()]
        pos_append = [[round(x, 4) for x in col] for col in mol_append.column_pos()]

        self.assertEqual(pos_gro, pos_pdb)
        self.assertEqual(pos_gro, pos_mol2)
        self.assertEqual(pos_gro, pos_atom)
        self.assertEqual([col+col for col in pos_gro], pos_concat)
        self.assertEqual([col+col for col in pos_gro], pos_append)

    def test_molecule_properties(self):
        mol = Molecule(inp="data/benzene.gro")

        self.assertEqual(mol.pos(0), [0.0935, 0.0000, 0.3143])
        self.assertEqual([[round(x, 4) for x in mol.bond(0, 1)[0]], round(mol.bond(0, 1)[1], 4)], [[0.1191, 0.0, 0.0687], 0.1375])
        self.assertEqual(mol.bond([1, 0, 0], [0, 0, 0]), [[-1, 0, 0], 1.0])
        self.assertEqual(mol.get_box(), [0.4252, 0.001, 0.491])
        self.assertEqual([round(x, 4) for x in mol.centroid()], [0.2126, 0.0, 0.2455])
        self.assertEqual([round(x, 4) for x in mol.com()], [0.2126, 0.0, 0.2455])

    def test_molecule_editing(self):
        mol = Molecule(inp="data/benzene.gro")

        mol.translate([0, 0.1, 0.2])
        self.assertEqual([round(x, 4) for x in mol.pos(3)], [0.3317, 0.1000, 0.3768])
        mol.rotate("x", 45)
        self.assertEqual([round(x, 4) for x in mol.pos(3)], [0.3317, -0.1957, 0.3371])
        mol.move(0, [1, 1, 1])
        self.assertEqual([round(x, 4) for x in mol.pos(3)], [1.2382, 1.0972, 0.9028])
        mol.zero()
        self.assertEqual([round(x, 4) for x in mol.pos(3)], [0.3317, 0.2222, 0.1250])
        mol.put(3, [0, 0, 0])
        self.assertEqual([round(x, 4) for x in mol.pos(3)], [0.0000, 0.0000, 0.0000])
        mol.part_move([0, 1], [1, 2, 3, 4], 0.5)
        self.assertEqual([round(x, 4) for x in mol.pos(3)], [0.3140, -0.1281, 0.1281])
        mol.part_rotate([0, 1], [1, 2, 3, 4], 45, 1)
        self.assertEqual([round(x, 4) for x in mol.pos(3)], [-0.1277, 0.0849, -0.3176])
        mol.part_angle([0, 1], [1, 2], [1, 2, 3, 4], 45, 1)
        self.assertEqual([round(x, 4) for x in mol.pos(3)], [-0.1360, -0.1084, -0.3068])

    def test_molecule_creation(self):
        mol = Molecule()

        mol.add("C", [0, 0.1, 0.2])
        mol.add("C", 0, r=0.1, theta=90)
        mol.add("C", 1, [0, 1], r=0.1, theta=90)
        mol.add("C", 2, [0, 2], r=0.1, theta=90, phi=45)
        self.assertEqual([round(x, 4) for x in mol.pos(3)], [0.0500, 0.0500, 0.0293])
        mol.delete(2)
        self.assertEqual([round(x, 4) for x in mol.pos(2)], [0.0500, 0.0500, 0.0293])
        mol.add("C", [0, 0.1, 0.2])
        self.assertEqual(mol.overlap(), {0: [3]})
        mol.switch_atom_order(0, 2)
        self.assertEqual([round(x, 4) for x in mol.pos(0)], [0.0500, 0.0500, 0.0293])
        mol.set_atom_type(0, "R")
        self.assertEqual(mol.get_atom_list()[0].get_atom_type(), "R")
        mol.set_atom_name(0, "RuX")
        self.assertEqual(mol.get_atom_list()[0].get_name(), "RuX")

    def test_molecule_set_get(self):
        mol = Molecule()

        mol.set_name("test_mol")
        mol.set_short("TMOL")
        mol.set_box([1, 1, 1])
        mol.set_charge(1.5)
        mol.set_masses([1, 2, 3])
        mol.set_mass(6)

        self.assertEqual(mol.get_name(), "test_mol")
        self.assertEqual(mol.get_short(), "TMOL")
        self.assertEqual(mol.get_box(), [1, 1, 1])
        self.assertEqual(mol.get_num(), 0)
        self.assertEqual(mol.get_charge(), 1.5)
        self.assertEqual(mol.get_masses(), [1, 2, 3])
        self.assertEqual(mol.get_mass(), 6)

    def test_molecule_representation(self):
        mol = Molecule()
        mol.add("H", [0.0, 0.1, 0.2], name="HO1")

        self.assertEqual(mol.__str__(), "  Name Type    x    y    z\n0  HO1    H  0.0  0.1  0.2")


    ##############
    # Essentails #
    ##############
    def test_essentials(self):
        self.assertEqual([round(x, 4) for x in Alkane(10).pos(5)], [0.0472, 0.1028, 0.7170])
        self.assertEqual([round(x, 4) for x in Alcohol(10).pos(5)], [0.0363, 0.1028, 0.7170])
        self.assertEqual([round(x, 4) for x in Ketone(10, 5).pos(5)], [0.0472, 0.1028, 0.7170])
        self.assertEqual([round(x, 4) for x in TMS(separation=30).pos(5)], [0.0273, 0.0472, 0.4525])


    #########
    # Store #
    #########
    def test_store(self):
        mol = Molecule(inp="data/benzene.gro")

        Store(mol, "output").gro()
        Store(mol, "output").pdb()
        Store(mol, "output").xyz()


    ###########
    # Pattern #
    ###########
    def test_pattern(self):
        beta_cristobalit = BetaCristobalit()

        # Pattern and output
        pattern = beta_cristobalit.pattern()
        pattern.set_name("beta_cristobalit_pattern")
        Store(pattern, "output").gro()
        self.assertEqual(pattern.get_num(), 36)

        # Generation and Orientation
        beta_cristobalit = BetaCristobalit()
        beta_cristobalit.generate([2, 2, 2], "x")
        self.assertEqual(beta_cristobalit.get_size(), [2.635, 1.827, 2.150])
        beta_cristobalit = BetaCristobalit()
        beta_cristobalit.generate([2, 2, 2], "y")
        self.assertEqual(beta_cristobalit.get_size(), [2.150, 2.635, 1.827])
        beta_cristobalit = BetaCristobalit()
        beta_cristobalit.generate([2, 2, 2], "z")
        self.assertEqual(beta_cristobalit.get_size(), [2.150, 1.827, 2.635])

        # Overlap and output
        block = beta_cristobalit.get_block()
        Store(block, "output").gro()
        self.assertEqual(block.get_num(), 576)
        self.assertEqual(block.overlap(), {})

        # Getter
        self.assertEqual(beta_cristobalit.get_repeat(), [0.506, 0.877, 1.240])
        self.assertEqual(beta_cristobalit.get_gap(), [0.126, 0.073, 0.155])


    ########
    # Cube #
    ########
    def test_cube(self):
        block = BetaCristobalit().generate([2, 2, 2], "z")
        block.set_name("cube")
        Store(block, "output").gro()
        cube = Cube(block, 0.4, True)

        # Splitting and filling
        self.assertEqual(len(cube.get_origin()), 120)
        self.assertEqual(cube.get_origin()[(1, 1, 1)], [0.4, 0.4, 0.4])
        self.assertEqual(cube.get_pointer()[(1, 1, 1)], [14, 46, 51, 52, 65])

        # Iterator
        self.assertEqual(cube._right((1, 1, 1)), (2, 1, 1))
        self.assertEqual(cube._left((1, 1, 1)),  (0, 1, 1))
        self.assertEqual(cube._top((1, 1, 1)),   (1, 2, 1))
        self.assertEqual(cube._bot((1, 1, 1)),   (1, 0, 1))
        self.assertEqual(cube._front((1, 1, 1)), (1, 1, 2))
        self.assertEqual(cube._back((1, 1, 1)),  (1, 1, 0))
        self.assertEqual(len(cube.neighbour((1, 1, 1))), 27)

        # Search
        self.assertEqual(cube.find_bond([(1, 1, 1)], ["Si", "O"], 0.155, 0.005), [[51, [14, 46, 52, 65]]])
        self.assertEqual(cube.find_bond([(1, 1, 1)], ["O", "Si"], 0.155, 0.005), [[14, [51, 13]], [46, [43, 51]], [52, [51, 49]], [65, [64, 51]]])
        self.assertEqual(cube.find_bond([(0, 0, 0)], ["Si", "O"], 0.155, 0.005), [[3, [4, 2, 174, 9]], [5, [306, 110, 4, 6]]])
        self.assertEqual(cube.find_bond([(0, 0, 0)], ["O", "Si"], 0.155, 0.005), [[4, [3, 5]], [6, [7, 5]]])

    def test_cube_parallel(self):
        self.skipTest("Parallel")

        cube = Cube(BetaCristobalit().generate([2, 2, 2], "z"), 0.4, True)

        self.assertEqual(len(cube.find_parallel(None, ["Si", "O"], 0.155, 0.005)), 192)
        self.assertEqual(len(cube.find_parallel(None, ["O", "Si"], 0.155, 0.005)), 384)


    ########
    # Pore #
    ########
    def test_pore(self):
        self.skipTest("Temporary")

        pore = Pore([2, 2, 2], "z")
        pore.generate(is_time=False)

        Store(pore.get_pore(), "output").gro()


if __name__ == '__main__':
    unittest.main(verbosity=2)
