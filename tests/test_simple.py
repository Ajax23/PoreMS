import os
import sys

import shutil
import unittest
import matplotlib.pyplot as plt

# # Install package
# if sys.platform == "win32":
#     os.system("pip install ../.")
# else:
#     os.system("pip install ../. &> /dev/null")
# print("Finished inistalling package...")

# Import package
import porems.utils as utils
import porems.database as db
import porems.geometry as geometry

from porems.atom import Atom
from porems.molecule import Molecule
from porems.essentials import *
from porems.store import Store
from porems.pattern import *
from porems.dice import Dice
from porems.matrix import Matrix
from porems.shape import *
from porems.pore import Pore


class UserModelCase(unittest.TestCase):
    #################
    # Remove Output #
    #################
    @classmethod
    def setUpClass(self):
        folder = 'output'
        if os.path.exists(folder):
            for filename in os.listdir(folder):
                file_path = os.path.join(folder, filename)
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
        else:
            os.makedirs(folder)


    #########
    # Utils #
    #########
    def test_utils(self):
        file_link = "output/test/test.txt"

        utils.mkdirp("output/test")

        with open(file_link, "w") as file_out:
            file_out.write("TEST")
        utils.copy(file_link, file_link+"t")
        utils.replace(file_link+"t", "TEST", "DOTA")
        with open(file_link+"t", "r") as file_in:
            for line in file_in:
                self.assertEqual(line, "DOTA\n")

        self.assertEqual(utils.column([[1, 1, 1], [2, 2, 2]]), [[1, 2], [1, 2], [1, 2]])

        utils.save([1, 1, 1], file_link)
        self.assertEqual(utils.load(file_link), [1, 1, 1])

        self.assertEqual(round(utils.toc(utils.tic(), is_print=True)), 0)


    ############
    # Geometry #
    ############
    def test_geometry(self):
        vec_a = [1, 1, 2]
        vec_b = [0, 3, 2]

        self.assertEqual(round(geometry.dot_product(vec_a, vec_b), 4), 7)
        self.assertEqual(round(geometry.length(vec_a), 4), 2.4495)
        self.assertEqual([round(x, 4) for x in geometry.vector(vec_a, vec_b)], [-1, 2, 0])
        self.assertIsNone(geometry.vector([0, 1], [0, 0, 0]))
        self.assertEqual([round(x, 4) for x in geometry.unit(vec_a)], [0.4082, 0.4082, 0.8165])
        self.assertEqual([round(x, 4) for x in geometry.cross_product(vec_a, vec_b)], [-4, -2, 3])
        self.assertEqual(round(geometry.angle(vec_a, vec_b), 4), 37.5714)
        self.assertEqual(round(geometry.angle_polar(vec_a), 4), 0.7854)
        self.assertEqual(round(geometry.angle_azi(vec_b), 4), 0.9828)
        self.assertEqual(round(geometry.angle_azi([0, 0, 0]), 4), 1.5708)
        self.assertEqual([round(x, 4) for x in geometry.main_axis(1)], [1, 0, 0])
        self.assertEqual([round(x, 4) for x in geometry.main_axis(2)], [0, 1, 0])
        self.assertEqual([round(x, 4) for x in geometry.main_axis(3)], [0, 0, 1])
        self.assertEqual([round(x, 4) for x in geometry.main_axis("x")], [1, 0, 0])
        self.assertEqual([round(x, 4) for x in geometry.main_axis("y")], [0, 1, 0])
        self.assertEqual([round(x, 4) for x in geometry.main_axis("z")], [0, 0, 1])
        self.assertEqual(geometry.main_axis("h"), "Wrong axis definition...")
        self.assertEqual(geometry.main_axis(100), "Wrong axis definition...")
        self.assertEqual(geometry.main_axis(0.1), "Wrong axis definition...")
        self.assertEqual([round(x, 4) for x in geometry.rotate(vec_a, "x", 90, True)], [1.0, -2.0, 1.0])
        self.assertIsNone(geometry.rotate(vec_a, [0, 1, 2, 3], 90, True))
        self.assertIsNone(geometry.rotate(vec_a, "h", 90, True))


    ############
    # Database #
    ############
    def test_database(self):
        self.assertEqual(db.get_mass("H"), 1.0079)
        self.assertIsNone(db.get_mass("DOTA"))


    ########
    # Atom #
    ########
    def test_atom(self):
        atom = Atom([0.0, 0.1, 0.2], "H", "HO1")

        self.assertEqual(atom.get_pos(), [0.0, 0.1, 0.2])
        self.assertEqual(atom.get_atom_type(), "H")
        self.assertEqual(atom.get_name(), "HO1")

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

        self.assertEqual(Molecule(inp="data/benzene.DOTA").get_atom_list(), None)

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
        mol.part_move([0, 1], [2, 3, 4], 0.5)
        mol.part_move([0, 1], 1, 0.5)
        self.assertEqual([round(x, 4) for x in mol.pos(3)], [0.3140, -0.1281, 0.1281])
        mol.part_rotate([0, 1], [2, 3, 4], 45, 1)
        mol.part_rotate([0, 1], 1, 45, 1)
        self.assertEqual([round(x, 4) for x in mol.pos(3)], [-0.1277, 0.0849, -0.3176])
        mol.part_angle([0, 1], [1, 2], [1, 2, 3, 4], 45, 1)
        self.assertEqual([round(x, 4) for x in mol.pos(3)], [-0.1360, -0.1084, -0.3068])
        mol.part_angle([0, 0, 1], [0, 1, 0], 1, 45, 1)
        self.assertEqual([round(x, 4) for x in mol.pos(3)], [-0.1360, -0.1084, -0.3068])

        self.assertIsNone(mol._vector(0.1, 0.1))
        self.assertIsNone(mol._vector([0, 0], [0, 0]))

        self.assertIsNone(mol.part_angle([0, 0, 1, 0], [0, 1, 0, 0], 1, 45, 1))
        self.assertIsNone(mol.part_angle([0, 0], [0, 1, 2], 1, 45, 1))

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
        self.assertEqual(mol.get_atom_type(0), "R")
        mol.set_atom_name(0, "RuX")
        self.assertEqual(mol.get_atom_list()[0].get_name(), "RuX")

    def test_molecule_set_get(self):
        mol = Molecule()

        mol.set_name("test_mol")
        mol.set_short("TMOL")
        mol.set_box([1, 1, 1])
        mol.set_charge(1.5)
        mol.set_masses([1, 2, 3])

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
        self.assertEqual([round(x, 4) for x in Alkane(10, "decane", "DEC").pos(5)], [0.0472, 0.1028, 0.7170])
        self.assertEqual([round(x, 4) for x in Alkane(1, "methane", "MET").pos(0)], [0.0514, 0.0890, 0.0363])
        self.assertEqual([round(x, 4) for x in Alcohol(10, "decanol", "DCOL").pos(5)], [0.0363, 0.1028, 0.7170])
        self.assertEqual([round(x, 4) for x in Alcohol(1, "methanol", "MEOL").pos(0)], [0.0715, 0.0890, 0.0363])
        self.assertEqual([round(x, 4) for x in Ketone(10, 5, "decanone", "DCON").pos(5)], [0.0472, 0.1028, 0.7170])
        self.assertEqual(Ketone(2, 0).get_name(), "ERROR")
        self.assertEqual([round(x, 4) for x in TMS(separation=30).pos(5)], [0.0273, 0.0472, 0.4525])
        self.assertEqual([round(x, 4) for x in TMS(is_si=False).pos(5)], [0.0273, 0.0472, 0.4976])


    #########
    # Store #
    #########
    def test_store(self):
        mol = Molecule(inp="data/benzene.gro")

        Store(mol, "output").job(True)
        Store(mol, "output").obj()
        Store(mol, "output").gro(use_atom_names=True)
        Store(mol, "output").pdb(use_atom_names=True)
        Store(mol, "output").xyz()
        Store(mol, "output").top()
        Store(mol, "output").grid()


    ###########
    # Pattern #
    ###########
    def test_pattern_beta_cristobalit(self):
        beta_cristobalit = BetaCristobalit()

        # Pattern and output
        pattern = beta_cristobalit.pattern()
        pattern.set_name("beta_cristobalit_pattern")
        self.assertEqual(pattern.get_num(), 36)
        Store(pattern, "output").gro()

        # Generation and Orientation
        beta_cristobalit = BetaCristobalit()
        beta_cristobalit.generate([2, 2, 2], "x")
        beta_cristobalit.get_block().set_name("beta_cristobalit_x")
        self.assertEqual(beta_cristobalit.get_size(), [2.480, 1.754, 2.024])
        self.assertEqual([round(x, 3) for x in beta_cristobalit.get_block().get_box()], [2.480, 1.754, 2.024])
        Store(beta_cristobalit.get_block(), "output").gro()

        beta_cristobalit = BetaCristobalit()
        beta_cristobalit.generate([2, 2, 2], "y")
        beta_cristobalit.get_block().set_name("beta_cristobalit_y")
        self.assertEqual(beta_cristobalit.get_size(), [2.024, 2.480, 1.754])
        self.assertEqual([round(x, 3) for x in beta_cristobalit.get_block().get_box()], [2.024, 2.480, 1.754])
        Store(beta_cristobalit.get_block(), "output").gro()

        beta_cristobalit = BetaCristobalit()
        beta_cristobalit.generate([2, 2, 2], "z")
        beta_cristobalit.get_block().set_name("beta_cristobalit_z")
        self.assertEqual(beta_cristobalit.get_size(), [2.024, 1.754, 2.480])
        self.assertEqual([round(x, 3) for x in beta_cristobalit.get_block().get_box()], [2.024, 1.754, 2.480])
        Store(beta_cristobalit.get_block(), "output").gro()

        # Exterior
        beta_cristobalit = BetaCristobalit()
        beta_cristobalit.generate([2, 2, 2], "x")
        beta_cristobalit.get_block().set_name("beta_cristobalit_ext_x")
        beta_cristobalit.exterior()
        self.assertEqual(beta_cristobalit.get_block().get_num(), 600)
        Store(beta_cristobalit.get_block(), "output").gro()

        beta_cristobalit = BetaCristobalit()
        beta_cristobalit.generate([2, 2, 2], "y")
        beta_cristobalit.get_block().set_name("beta_cristobalit_ext_y")
        self.assertIsNone(beta_cristobalit.exterior())

        beta_cristobalit = BetaCristobalit()
        beta_cristobalit.generate([2, 2, 2], "z")
        beta_cristobalit.get_block().set_name("beta_cristobalit_ext_z")
        beta_cristobalit.exterior()
        self.assertEqual(beta_cristobalit.get_block().get_num(), 592)
        Store(beta_cristobalit.get_block(), "output").gro()

        # Misc
        beta_cristobalit = BetaCristobalit()
        beta_cristobalit.generate([2, 2, 2], "z")
        beta_cristobalit.get_block().set_name("test_mol")

        # Overlap and output
        self.assertEqual(beta_cristobalit.get_block().get_num(), 576)
        self.assertEqual(beta_cristobalit.get_block().overlap(), {})

        # Getter
        self.assertEqual(beta_cristobalit.get_repeat(), [0.506, 0.877, 1.240])
        self.assertEqual(beta_cristobalit.get_gap(), [0.126, 0.073, 0.155])
        self.assertEqual(beta_cristobalit.get_orient(), "z")
        self.assertEqual(beta_cristobalit.get_block().get_name(), "test_mol")


    ########
    # Dice #
    ########
    def test_dice(self):
        block = BetaCristobalit().generate([2, 2, 2], "z")
        block.set_name("dice")
        Store(block, "output").gro()
        dice = Dice(block, 0.4, True)

        # Splitting and filling
        self.assertEqual(len(dice.get_origin()), 120)
        self.assertEqual(dice.get_origin()[(1, 1, 1)], [0.4, 0.4, 0.4])
        self.assertEqual(dice.get_pointer()[(1, 1, 1)], [14, 46, 51, 52, 65])

        # Iterator
        self.assertEqual(dice._right((1, 1, 1)), (2, 1, 1))
        self.assertEqual(dice._left((1, 1, 1)),  (0, 1, 1))
        self.assertEqual(dice._top((1, 1, 1)),   (1, 2, 1))
        self.assertEqual(dice._bot((1, 1, 1)),   (1, 0, 1))
        self.assertEqual(dice._front((1, 1, 1)), (1, 1, 2))
        self.assertEqual(dice._back((1, 1, 1)),  (1, 1, 0))
        self.assertEqual(len(dice.neighbor((1, 1, 1))), 27)
        self.assertEqual(len(dice.neighbor((1, 1, 1), False)), 26)

        # Search
        self.assertEqual(dice.find_bond([(1, 1, 1)], ["Si", "O"], 0.155, 0.005), [[51, [14, 46, 52, 65]]])
        self.assertEqual(dice.find_bond([(1, 1, 1)], ["O", "Si"], 0.155, 0.005), [[14, [51, 13]], [46, [43, 51]], [52, [51, 49]], [65, [64, 51]]])
        self.assertEqual(dice.find_bond([(0, 0, 0)], ["Si", "O"], 0.155, 0.005), [[3, [4, 2, 174, 9]], [5, [306, 110, 4, 6]]])
        self.assertEqual(dice.find_bond([(0, 0, 0)], ["O", "Si"], 0.155, 0.005), [[4, [3, 5]], [6, [7, 5]]])

        # Parallel search
        self.assertEqual(len(dice.find_parallel(None, ["Si", "O"], 0.155, 0.005)), 192)
        self.assertEqual(len(dice.find_parallel(None, ["O", "Si"], 0.155, 0.005)), 384)

        # Setter Getter
        dice.set_pbc(True)
        self.assertEqual(dice.get_count(), [5, 4, 6])
        self.assertEqual(dice.get_size(), 0.4)
        self.assertEqual(dice.get_mol().get_name(), "dice")


    ##########
    # Matrix #
    ##########
    def test_matrix(self):
        orient = "z"
        block = BetaCristobalit().generate([1, 1, 1], orient)
        block.set_name("matrix")
        Store(block, "output").gro()
        dice = Dice(block, 0.2, True)
        bonds = dice.find_bond(None, ["Si", "O"], 0.155, 10e-2)

        matrix = Matrix(bonds)
        connect = matrix.get_matrix()
        matrix.split(0, 17)
        self.assertEqual(connect[0], [30, 8, 1])
        self.assertEqual(connect[17], [19])
        matrix.strip(0)
        self.assertEqual(connect[0], [])
        self.assertEqual(connect[1], [43])
        self.assertEqual(connect[8], [7])
        self.assertEqual(connect[30], [3])
        self.assertEqual(matrix.bound(0), [0])
        self.assertEqual(matrix.bound(1, "lt"), [0])
        self.assertEqual(matrix.bound(4, "gt"), [])
        self.assertIsNone(matrix.bound(4, "test"))


    #########
    # Shape #
    #########
    def test_shape_cylinder(self):
        self.skipTest("Temporary")

        block = BetaCristobalit().generate([6, 6, 6], "z")
        block.set_name("shape_cylinder")
        dice = Dice(block, 0.4, True)
        matrix = Matrix(dice.find_parallel(None, ["Si", "O"], 0.155, 10e-2))
        centroid = block.centroid()
        central = geometry.unit(geometry.rotate([0, 0, 1], [1, 0, 0], 45, True))

        cylinder = Cylinder({"centroid": centroid, "central": central, "length": 3, "diameter": 4})

        # Properties
        self.assertEqual(round(cylinder.volume(), 4), 37.6991)
        self.assertEqual(round(cylinder.surface(), 4), 37.6991)

        # Test vector
        vec = [3.6716, 4.4441, 0.2840]

        # Surface
        self.assertEqual([round(x[0][20], 4) for x in cylinder.surf(num=100)], vec)
        self.assertEqual([round(x[0][20], 4) for x in cylinder.rim(0, num=100)], vec)

        # Normal
        self.assertEqual([round(x, 4) for x in cylinder.convert([0, 0, 0], False)], [3.0777, 3.0937, 1.6344])
        self.assertEqual([round(x, 4) for x in cylinder.normal(vec)], [0.5939, 2.9704, 0.0000])

        # Positioning
        del_list = [atom_id for atom_id, atom in enumerate(block.get_atom_list()) if cylinder.is_in(atom.get_pos())]
        matrix.strip(del_list)
        block.delete(matrix.bound(0))
        self.assertEqual(block.get_num(), 12650)

        # Store molecule
        Store(block, "output").gro()

        # Plot surface
        cylinder.plot(vec=[3.17290646, 4.50630614, 0.22183271])
        # plt.show()

    def test_shape_sphere(self):
        self.skipTest("Temporary")

        block = BetaCristobalit().generate([6, 6, 6], "z")
        block.set_name("shape_sphere")
        dice = Dice(block, 0.4, True)
        matrix = Matrix(dice.find_parallel(None, ["Si", "O"], 0.155, 10e-2))
        centroid = block.centroid()
        central = geometry.unit(geometry.rotate([0, 0, 1], [1, 0, 0], 0, True))

        sphere = Sphere({"centroid": centroid, "central": central, "diameter": 4})

        # Properties
        self.assertEqual(round(sphere.volume(), 4), 33.5103)
        self.assertEqual(round(sphere.surface(), 4), 50.2655)

        # Surface
        self.assertEqual([round(x[0][20], 4) for x in sphere.surf(num=100)], [4.2636, 3.0937, 4.745])
        self.assertEqual([round(x[0][20], 4) for x in sphere.rim(0, num=100)], [4.9875, 3.0937, 3.7283])

        # Normal
        self.assertEqual([round(x, 4) for x in sphere.convert([0, 0, 0], False)], [3.0777, 3.0937, 3.1344])
        self.assertEqual([round(x, 4) for x in sphere.normal([4.2636, 3.0937, 4.745])], [1.4063, 0.0000, 1.9099])

        # Positioning
        del_list = [atom_id for atom_id, atom in enumerate(block.get_atom_list()) if sphere.is_in(atom.get_pos())]
        matrix.strip(del_list)
        block.delete(matrix.bound(0))
        self.assertEqual(block.get_num(), 12934)

        # Store molecule
        Store(block, "output").gro()

        # Plot surface
        sphere.plot(inp=3.14, vec=[1.08001048, 3.09687610, 1.72960828])
        # plt.show()

    def test_shape_capsule(self):
        self.skipTest("Temporary")

        pattern = BetaCristobalit()
        block = pattern.generate([6, 6, 12], "z")
        block.set_name("shape_capsule")
        dice = Dice(block, 0.4, True)
        matrix = Matrix(dice.find_parallel(None, ["Si", "O"], 0.155, 10e-2))
        central = geometry.unit(geometry.rotate([0, 0, 1], [1, 0, 0], 0, True))
        centroid = block.centroid()
        centroid_cyl_l = centroid[:2]+[0]
        centroid_cyl_r = centroid[:2]+[pattern.get_size()[2]]
        centroid_sph_l = centroid[:2]+[3]
        centroid_sph_r = centroid[:2]+[pattern.get_size()[2]-3]

        cylinder_l = Cylinder({"centroid": centroid_cyl_l, "central": central, "length": 6, "diameter": 4})
        cylinder_r = Cylinder({"centroid": centroid_cyl_r, "central": central, "length": 6, "diameter": 4})
        sphere_l = Sphere({"centroid": centroid_sph_l, "central": central, "diameter": 4})
        sphere_r = Sphere({"centroid": centroid_sph_r, "central": central, "diameter": 4})

        del_list = []
        del_list.extend([atom_id for atom_id, atom in enumerate(block.get_atom_list()) if cylinder_l.is_in(atom.get_pos())])
        del_list.extend([atom_id for atom_id, atom in enumerate(block.get_atom_list()) if cylinder_r.is_in(atom.get_pos())])
        del_list.extend([atom_id for atom_id, atom in enumerate(block.get_atom_list()) if sphere_l.is_in(atom.get_pos())])
        del_list.extend([atom_id for atom_id, atom in enumerate(block.get_atom_list()) if sphere_r.is_in(atom.get_pos())])
        matrix.strip(del_list)

        block.delete(matrix.bound(0))
        Store(block, "output").gro()


    ########
    # Pore #
    ########
    def test_pore(self):
        orient = "z"
        pattern = BetaCristobalit()
        pattern.generate([6, 6, 6], orient)
        pattern.exterior()

        block = pattern.get_block()
        block.set_name("pore")

        dice = Dice(block, 0.4, True)
        matrix = Matrix(dice.find_parallel(None, ["Si", "O"], 0.155, 10e-2))
        oxygen_out = matrix.bound(1)

        centroid = block.centroid()
        central = geometry.unit(geometry.rotate([0, 0, 1], [1, 0, 0], 0, True))
        cylinder = Cylinder({"centroid": centroid, "central": central, "length": 6, "diameter": 4})
        del_list = [atom_id for atom_id, atom in enumerate(block.get_atom_list()) if cylinder.is_in(atom.get_pos())]
        matrix.strip(del_list)

        # Create pore object
        pore = Pore(block, matrix)
        pore.prepare()
        self.assertEqual(len(matrix.bound(1)), 710)
        pore.sites(oxygen_out)

        sites = pore.get_sites()

        # Output
        block.delete(matrix.bound(0))
        Store(block, "output").gro()


if __name__ == '__main__':
    unittest.main(verbosity=2)
