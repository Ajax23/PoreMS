import os
import sys

import shutil
import unittest
import matplotlib.pyplot as plt

import porems as pms


class UserModelCase(unittest.TestCase):
    #################
    # Remove Output #
    #################
    @classmethod
    def setUpClass(self):
        if os.path.isdir("tests"):
            os.chdir("tests")

        folder = 'output'
        pms.utils.mkdirp(folder)
        pms.utils.mkdirp(folder+"/temp")
        open(folder+"/temp.txt", 'a').close()

        for filename in os.listdir(folder):
            file_path = os.path.join(folder, filename)
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)


    #########
    # Utils #
    #########
    def test_utils(self):
        file_link = "output/test/test.txt"

        pms.utils.mkdirp("output/test")

        with open(file_link, "w") as file_out:
            file_out.write("TEST")
        pms.utils.copy(file_link, file_link+"t")
        pms.utils.replace(file_link+"t", "TEST", "DOTA")
        with open(file_link+"t", "r") as file_in:
            for line in file_in:
                self.assertEqual(line, "DOTA\n")

        self.assertEqual(pms.utils.column([[1, 1, 1], [2, 2, 2]]), [[1, 2], [1, 2], [1, 2]])

        pms.utils.save([1, 1, 1], file_link)
        self.assertEqual(pms.utils.load(file_link), [1, 1, 1])

        self.assertEqual(round(pms.utils.mumol_m2_to_mols(3, 100), 4), 180.66)
        self.assertEqual(round(pms.utils.mols_to_mumol_m2(180, 100), 4), 2.989)
        self.assertEqual(round(pms.utils.mmol_g_to_mumol_m2(0.072, 512), 2), 0.14)
        self.assertEqual(round(pms.utils.mmol_l_to_mols(30, 1000), 4), 18.066)
        self.assertEqual(round(pms.utils.mols_to_mmol_l(18, 1000), 4), 29.8904)

        print()
        pms.utils.toc(pms.utils.tic(), message="Test", is_print=True)
        self.assertEqual(round(pms.utils.toc(pms.utils.tic(), is_print=True)), 0)


    ############
    # Geometry #
    ############
    def test_geometry(self):
        vec_a = [1, 1, 2]
        vec_b = [0, 3, 2]

        print()

        self.assertEqual(round(pms.geom.dot_product(vec_a, vec_b), 4), 7)
        self.assertEqual(round(pms.geom.length(vec_a), 4), 2.4495)
        self.assertEqual([round(x, 4) for x in pms.geom.vector(vec_a, vec_b)], [-1, 2, 0])
        self.assertIsNone(pms.geom.vector([0, 1], [0, 0, 0]))
        self.assertEqual([round(x, 4) for x in pms.geom.unit(vec_a)], [0.4082, 0.4082, 0.8165])
        self.assertEqual([round(x, 4) for x in pms.geom.cross_product(vec_a, vec_b)], [-4, -2, 3])
        self.assertEqual(round(pms.geom.angle(vec_a, vec_b), 4), 37.5714)
        self.assertEqual(round(pms.geom.angle_polar(vec_a), 4), 0.7854)
        self.assertEqual(round(pms.geom.angle_azi(vec_b), 4), 0.9828)
        self.assertEqual(round(pms.geom.angle_azi([0, 0, 0]), 4), 1.5708)
        self.assertEqual([round(x, 4) for x in pms.geom.main_axis(1)], [1, 0, 0])
        self.assertEqual([round(x, 4) for x in pms.geom.main_axis(2)], [0, 1, 0])
        self.assertEqual([round(x, 4) for x in pms.geom.main_axis(3)], [0, 0, 1])
        self.assertEqual([round(x, 4) for x in pms.geom.main_axis("x")], [1, 0, 0])
        self.assertEqual([round(x, 4) for x in pms.geom.main_axis("y")], [0, 1, 0])
        self.assertEqual([round(x, 4) for x in pms.geom.main_axis("z")], [0, 0, 1])
        self.assertEqual(pms.geom.main_axis("h"), "Wrong axis definition...")
        self.assertEqual(pms.geom.main_axis(100), "Wrong axis definition...")
        self.assertEqual(pms.geom.main_axis(0.1), "Wrong axis definition...")
        self.assertEqual([round(x, 4) for x in pms.geom.rotate(vec_a, "x", 90, True)], [1.0, -2.0, 1.0])
        self.assertIsNone(pms.geom.rotate(vec_a, [0, 1, 2, 3], 90, True))
        self.assertIsNone(pms.geom.rotate(vec_a, "h", 90, True))


    ############
    # Database #
    ############
    def test_database(self):
        print()
        self.assertEqual(pms.db.get_mass("H"), 1.0079)
        self.assertIsNone(pms.db.get_mass("DOTA"))


    ########
    # Atom #
    ########
    def test_atom(self):
        atom = pms.Atom([0, 0, 0], "O", "O", 5)

        atom.set_pos([0.0, 0.1, 0.2])
        atom.set_atom_type("H")
        atom.set_name("HO1")
        atom.set_residue(0)

        self.assertEqual(atom.get_pos(), [0.0, 0.1, 0.2])
        self.assertEqual(atom.get_atom_type(), "H")
        self.assertEqual(atom.get_name(), "HO1")
        self.assertEqual(atom.get_residue(), 0)

        self.assertEqual(atom.__str__(), "   Residue Name Type    x    y    z\n0        0  HO1    H  0.0  0.1  0.2")


    ############
    # Molecule #
    ############
    def test_molecule_loading(self):
        mol_gro = pms.Molecule(inp="data/benzene.gro")
        mol_pdb = pms.Molecule(inp="data/benzene.pdb")
        mol_mol2 = pms.Molecule(inp="data/benzene.mol2")

        mol_atom = pms.Molecule(inp=mol_mol2.get_atom_list())
        mol_concat = pms.Molecule(inp=[mol_gro, mol_pdb])

        mol_append = pms.Molecule(inp="data/benzene.gro")
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

        print()
        self.assertEqual(pms.Molecule(inp="data/benzene.DOTA").get_atom_list(), None)

    def test_molecule_properties(self):
        mol = pms.Molecule(inp="data/benzene.gro")

        self.assertEqual(mol.pos(0), [0.0935, 0.0000, 0.3143])
        self.assertEqual([round(x, 4) for x in mol.bond(0, 1)], [0.1191, 0.0, 0.0687])
        self.assertEqual(mol.bond([1, 0, 0], [0, 0, 0]), [-1, 0, 0])
        self.assertEqual(mol.get_box(), [0.4252, 0.001, 0.491])
        self.assertEqual([round(x, 4) for x in mol.centroid()], [0.2126, 0.0, 0.2455])
        self.assertEqual([round(x, 4) for x in mol.com()], [0.2126, 0.0, 0.2455])

    def test_molecule_editing(self):
        mol = pms.Molecule(inp="data/benzene.gro")

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

        print()

        self.assertIsNone(mol._vector(0.1, 0.1))
        self.assertIsNone(mol._vector([0, 0], [0, 0]))

        self.assertIsNone(mol.part_angle([0, 0, 1, 0], [0, 1, 0, 0], 1, 45, 1))
        self.assertIsNone(mol.part_angle([0, 0], [0, 1, 2], 1, 45, 1))

    def test_molecule_creation(self):
        mol = pms.Molecule()

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
        mol.set_atom_residue(0, 1)
        self.assertEqual(mol.get_atom_list()[0].get_residue(), 1)

    def test_molecule_set_get(self):
        mol = pms.Molecule()

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
        mol = pms.Molecule()
        mol.add("H", [0.0, 0.1, 0.2], name="HO1")

        self.assertEqual(mol.__str__(), "   Residue Name Type    x    y    z\n0        0  HO1    H  0.0  0.1  0.2")


    ###########
    # Generic #
    ###########
    def test_generic(self):
        self.assertEqual([round(x, 4) for x in pms.gen.alkane(10, "decane", "DEC").pos(5)], [0.0472, 0.1028, 0.7170])
        self.assertEqual([round(x, 4) for x in pms.gen.alkane(1, "methane", "MET").pos(0)], [0.0514, 0.0890, 0.0363])
        self.assertEqual([round(x, 4) for x in pms.gen.alcohol(10, "decanol", "DCOL").pos(5)], [0.0363, 0.1028, 0.7170])
        self.assertEqual([round(x, 4) for x in pms.gen.alcohol(1, "methanol", "MEOL").pos(0)], [0.0715, 0.0890, 0.0363])
        self.assertEqual([round(x, 4) for x in pms.gen.ketone(10, 5, "decanone", "DCON").pos(5)], [0.0472, 0.1028, 0.7170])
        self.assertEqual([round(x, 4) for x in  pms.gen.tms(separation=30).pos(5)], [0.0273, 0.0472, 0.4525])
        self.assertEqual([round(x, 4) for x in  pms.gen.tms(is_si=False).pos(5)], [0.0273, 0.0472, 0.4976])
        self.assertEqual([round(x, 4) for x in  pms.gen.silanol().pos(0)], [0.000, 0.000, 0.000])

        print()
        self.assertIsNone(pms.gen.ketone(2, 0))


    #########
    # Store #
    #########
    def test_store(self):
        mol = pms.Molecule(inp="data/benzene.gro")

        mol.set_atom_residue(1, 1)

        pms.Store(mol, "output").job("store_job", "store_master.job")
        pms.Store(mol, "output").obj("store_obj.obj")
        pms.Store(mol, "output").gro("store_gro.gro", True)
        pms.Store(mol, "output").pdb("store_pdb.pdb", True)
        pms.Store(mol, "output").xyz("store_xyz.xyz")
        pms.Store(mol, "output").lmp("store_lmp.lmp")
        pms.Store(mol, "output").grid("store_grid.itp")

        print()
        pms.Store({})
        self.assertIsNone(pms.Store(mol).top())


    ###########
    # Pattern #
    ###########
    def test_pattern_beta_cristobalit(self):
        # Initialize
        beta_cristobalit = pms.BetaCristobalit()

        # Pattern and output
        pattern = beta_cristobalit.pattern()
        pattern.set_name("pattern_beta_cbt_minimal")
        self.assertEqual(pattern.get_num(), 36)
        pms.Store(pattern, "output").gro()

        # Generation and Orientation
        beta_cristobalit = pms.BetaCristobalit()
        beta_cristobalit.generate([2, 2, 2], "x")
        beta_cristobalit.get_block().set_name("pattern_beta_cbt_x")
        self.assertEqual(beta_cristobalit.get_size(), [2.480, 1.754, 2.024])
        self.assertEqual([round(x, 3) for x in beta_cristobalit.get_block().get_box()], [2.480, 1.754, 2.024])
        pms.Store(beta_cristobalit.get_block(), "output").gro()

        beta_cristobalit = pms.BetaCristobalit()
        beta_cristobalit.generate([2, 2, 2], "y")
        beta_cristobalit.get_block().set_name("pattern_beta_cbt_y")
        self.assertEqual(beta_cristobalit.get_size(), [2.024, 2.480, 1.754])
        self.assertEqual([round(x, 3) for x in beta_cristobalit.get_block().get_box()], [2.024, 2.480, 1.754])
        pms.Store(beta_cristobalit.get_block(), "output").gro()

        beta_cristobalit = pms.BetaCristobalit()
        beta_cristobalit.generate([2, 2, 2], "z")
        beta_cristobalit.get_block().set_name("pattern_beta_cbt_z")
        self.assertEqual(beta_cristobalit.get_size(), [2.024, 1.754, 2.480])
        self.assertEqual([round(x, 3) for x in beta_cristobalit.get_block().get_box()], [2.024, 1.754, 2.480])
        pms.Store(beta_cristobalit.get_block(), "output").gro()
        pms.Store(beta_cristobalit.get_block(), "output").lmp()

        # Misc
        beta_cristobalit = pms.BetaCristobalit()
        beta_cristobalit.generate([2, 2, 2], "z")
        beta_cristobalit.get_block().set_name("DOTA")

        # Overlap and output
        self.assertEqual(beta_cristobalit.get_block().get_num(), 576)
        self.assertEqual(beta_cristobalit.get_block().overlap(), {})

        # Getter
        self.assertEqual(beta_cristobalit.get_repeat(), [0.506, 0.877, 1.240])
        self.assertEqual(beta_cristobalit.get_gap(), [0.126, 0.073, 0.155])
        self.assertEqual(beta_cristobalit.get_orient(), "z")
        self.assertEqual(beta_cristobalit.get_block().get_name(), "DOTA")


    ########
    # Dice #
    ########
    def test_dice(self):
        block = pms.BetaCristobalit().generate([2, 2, 2], "z")
        block.set_name("dice")
        pms.Store(block, "output").gro()
        dice = pms.Dice(block, 0.4, True)

        # Splitting and filling
        self.assertEqual(len(dice.get_origin()), 120)
        self.assertEqual(dice.get_origin()[(1, 1, 1)], [0.4, 0.4, 0.4])
        self.assertEqual(dice.get_pointer()[(1, 1, 1)], [14, 51, 52, 64, 65, 67])

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
        self.assertEqual(dice.find_bond([(1, 1, 1)], ["Si", "O"], [0.155-0.005, 0.155+0.005]), [[51, [46, 14, 52, 65]], [64, [26, 63, 65, 67]]])
        self.assertEqual(dice.find_bond([(1, 1, 1)], ["O", "Si"], [0.155-0.005, 0.155+0.005]), [[14, [51, 13]], [52, [51, 49]], [65, [51, 64]], [67, [64, 69]]])
        self.assertEqual(dice.find_bond([(0, 0, 0)], ["Si", "O"], [0.155-0.005, 0.155+0.005]), [[3, [4, 9, 2, 174]], [5, [306, 110, 4, 6]]])
        self.assertEqual(dice.find_bond([(0, 0, 0)], ["O", "Si"], [0.155-0.005, 0.155+0.005]), [[4, [3, 5]], [6, [7, 5]], [9, [3, 11]]])

        # Parallel search
        self.assertEqual(len(dice.find_parallel(None, ["Si", "O"], [0.155-0.005, 0.155+0.005])), 192)
        self.assertEqual(len(dice.find_parallel(None, ["O", "Si"], [0.155-0.005, 0.155+0.005])), 384)

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
        block = pms.BetaCristobalit().generate([1, 1, 1], orient)
        block.set_name("matrix")
        pms.Store(block, "output").gro()
        dice = pms.Dice(block, 0.2, True)
        bonds = dice.find_bond(None, ["Si", "O"], [0.155-1e-2, 0.155+1e-2])

        matrix = pms.Matrix(bonds)
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
        matrix.add(0, 17)
        self.assertEqual(connect[0], [17])
        self.assertEqual(connect[17], [19, 0])

        print()
        self.assertIsNone(matrix.bound(4, "test"))


    #########
    # Shape #
    #########
    def test_shape_cylinder(self):
        # self.skipTest("Temporary")

        block = pms.BetaCristobalit().generate([6, 6, 6], "z")
        block.set_name("shape_cylinder")
        dice = pms.Dice(block, 0.4, True)
        matrix = pms.Matrix(dice.find_parallel(None, ["Si", "O"], [0.155-1e-2, 0.155+1e-2]))
        centroid = block.centroid()
        central = pms.geom.unit(pms.geom.rotate([0, 0, 1], [1, 0, 0], 45, True))

        cylinder = pms.Cylinder({"centroid": centroid, "central": central, "length": 3, "diameter": 4})

        # Properties
        self.assertEqual(round(cylinder.volume(), 4), 37.6991)
        self.assertEqual(round(cylinder.surface(), 4), 37.6991)

        # Test vector
        vec = [3.6086, 4.4076, 0.2065]

        # Surface
        self.assertEqual([round(x[0][20], 4) for x in cylinder.surf(num=100)], vec)
        self.assertEqual([round(x[0][20], 4) for x in cylinder.rim(0, num=100)], vec)

        # Normal
        self.assertEqual([round(x, 4) for x in cylinder.convert([0, 0, 0], False)], [3.0147, 3.0572, 1.5569])
        self.assertEqual([round(x, 4) for x in cylinder.normal(vec)], [0.5939, 2.9704, 0.0000])

        # Positioning
        del_list = [atom_id for atom_id, atom in enumerate(block.get_atom_list()) if cylinder.is_in(atom.get_pos())]
        matrix.strip(del_list)
        block.delete(matrix.bound(0))
        self.assertEqual(block.get_num(), 12650)

        # Store molecule
        pms.Store(block, "output").gro()

        # Plot surface
        plt.figure()
        cylinder.plot(vec=[3.17290646, 4.50630614, 0.22183271])
        # plt.show()

    def test_shape_sphere(self):
        # self.skipTest("Temporary")

        block = pms.BetaCristobalit().generate([6, 6, 6], "z")
        block.set_name("shape_sphere")
        dice = pms.Dice(block, 0.4, True)
        matrix = pms.Matrix(dice.find_parallel(None, ["Si", "O"], [0.155-1e-2, 0.155+1e-2]))
        centroid = block.centroid()
        central = pms.geom.unit(pms.geom.rotate([0, 0, 1], [1, 0, 0], 0, True))

        sphere = pms.Sphere({"centroid": centroid, "central": central, "diameter": 4})

        # Properties
        self.assertEqual(round(sphere.volume(), 4), 33.5103)
        self.assertEqual(round(sphere.surface(), 4), 50.2655)

        # Surface
        self.assertEqual([round(x[0][20], 4) for x in sphere.surf(num=100)], [4.2006, 3.0572, 4.6675])
        self.assertEqual([round(x[0][20], 4) for x in sphere.rim(0, num=100)], [4.9245, 3.0572, 3.6508])

        # Normal
        self.assertEqual([round(x, 4) for x in sphere.convert([0, 0, 0], False)], [3.0147, 3.0572, 3.0569])
        self.assertEqual([round(x, 4) for x in sphere.normal([4.2006, 3.0572, 4.6675])], [1.4063, 0.0000, 1.9099])

        # Positioning
        del_list = [atom_id for atom_id, atom in enumerate(block.get_atom_list()) if sphere.is_in(atom.get_pos())]
        matrix.strip(del_list)
        block.delete(matrix.bound(0))
        self.assertEqual(block.get_num(), 12934)

        # Store molecule
        pms.Store(block, "output").gro()

        # Plot surface
        sphere.plot(inp=3.14, vec=[1.08001048, 3.09687610, 1.72960828])
        # plt.show()

    def test_shape_cuboid(self):
        # self.skipTest("Temporary")

        block = pms.BetaCristobalit().generate([6, 6, 6], "z")
        block.set_name("shape_cuboid")
        dice = pms.Dice(block, 0.4, True)
        matrix = pms.Matrix(dice.find_parallel(None, ["Si", "O"], [0.155-1e-2, 0.155+1e-2]))
        centroid = block.centroid()
        central = pms.geom.unit(pms.geom.rotate([0, 0, 1], [1, 0, 0], 0, True))

        cuboid = pms.Cuboid({"centroid": centroid, "central": central, "length": 10, "width": 6, "height": 4})

        # Properties
        self.assertEqual(round(cuboid.volume(), 4), 240)
        self.assertEqual(round(cuboid.surface(), 4), 248)

        # Normal
        self.assertEqual([round(x, 4) for x in cuboid.convert([0, 0, 0], False)], [0.0147, 1.0572, -1.9431])
        self.assertEqual([round(x, 4) for x in cuboid.normal([4.2636, 3.0937, 4.745])], [0, 1, 0])

        # Positioning
        del_list = [atom_id for atom_id, atom in enumerate(block.get_atom_list()) if cuboid.is_in(atom.get_pos())]
        matrix.strip(del_list)
        block.delete(matrix.bound(0))
        self.assertEqual(block.get_num(), 5160)

        # Store molecule
        pms.Store(block, "output").gro()

        # Plot surface
        cuboid.plot()
        # plt.show()


    ########
    # Pore #
    ########
    def test_pore(self):
        # self.skipTest("Temporary")

        # No exterior surface
        orient = "z"
        pattern = pms.BetaCristobalit()
        pattern.generate([6, 6, 6], orient)

        block = pattern.get_block()
        block.set_name("pore_cylinder_block")

        dice = pms.Dice(block, 0.4, True)
        bond_list = dice.find_parallel(None, ["Si", "O"], [0.155-1e-2, 0.155+1e-2])
        matrix = pms.Matrix(bond_list)

        pore = pms.Pore(block, matrix)

        centroid = block.centroid()
        central = pms.geom.unit(pms.geom.rotate([0, 0, 1], [1, 0, 0], 0, True))
        cylinder = pms.Cylinder({"centroid": centroid, "central": central, "length": 6, "diameter": 4})
        del_list = [atom_id for atom_id, atom in enumerate(block.get_atom_list()) if cylinder.is_in(atom.get_pos())]
        matrix.strip(del_list)

        pore.prepare()
        pore.sites()
        self.assertEqual(len(pore.get_sites()), 455)

        block.delete(matrix.bound(0))
        pms.Store(block, "output").gro("pore_no_ex.gro")

        # With exterior surface
        orient = "z"
        pattern = pms.BetaCristobalit()
        pattern.generate([6, 6, 6], orient)

        block = pattern.get_block()
        block.set_name("pore_cylinder_block")

        dice = pms.Dice(block, 0.4, True)
        bond_list = dice.find_parallel(None, ["Si", "O"], [0.155-1e-2, 0.155+1e-2])
        matrix = pms.Matrix(bond_list)

        pore = pms.Pore(block, matrix)
        pore.exterior(pattern.get_gap())

        centroid = block.centroid()
        central = pms.geom.unit(pms.geom.rotate([0, 0, 1], [1, 0, 0], 0, True))
        cylinder = pms.Cylinder({"centroid": centroid, "central": central, "length": 6, "diameter": 4})
        del_list = [atom_id for atom_id, atom in enumerate(block.get_atom_list()) if cylinder.is_in(atom.get_pos())]
        matrix.strip(del_list)

        pore.prepare()
        pore.amorph()
        self.assertEqual(len(matrix.bound(1)), 710)
        pore.sites()
        site_list = pore.get_sites()
        site_in = [site_key for site_key, site_val in site_list.items() if site_val["type"]=="in"]
        site_ex = [site_key for site_key, site_val in site_list.items() if site_val["type"]=="ex"]
        self.assertEqual(len(site_in), 432)
        self.assertEqual(len(site_ex), 201)

        si_pos_in = [block.pos(site_key) for site_key, site_val in site_list.items() if site_val["type"]=="in"]
        si_pos_ex = [block.pos(site_key) for site_key, site_val in site_list.items() if site_val["type"]=="ex"]

        if si_pos_in:
            temp_mol = pms.Molecule()
            for pos in si_pos_in:
                temp_mol.add("Si", pos)
            size = temp_mol.get_box()[2]
            pms.Store(temp_mol).gro("output/pore_cylinder_si_in.gro")

        if si_pos_ex:
            temp_mol = pms.Molecule()
            for pos in si_pos_ex:
                temp_mol.add("Si", pos)
            size = temp_mol.get_box()[2]
            pms.Store(temp_mol).gro("output/pore_cylinder_si_ex.gro")

        # Objectify grid
        non_grid = matrix.bound(1)+list(site_list.keys())
        bonded = matrix.bound(0, "gt")
        grid_atoms = [atom for atom in bonded if not atom in non_grid]
        mol_obj = pore.objectify(grid_atoms)
        self.assertEqual(len(mol_obj), 8279)
        pms.Store(pms.Molecule(name="pore_cylinder_grid", inp=mol_obj), "output").gro(use_atom_names=True)

        # Attachment
        mol = pms.gen.tms()

        def normal(pos):
            return [0, 0, -1] if pos[2] < centroid[2] else [0, 0, 1]

        ## Siloxane
        mols_siloxane = pore.siloxane(site_in, 100, cylinder.normal)
        site_in = [site_key for site_key, site_val in site_list.items() if site_val["type"]=="in"]

        ## Normal
        mols_in = pore.attach(mol, 0, [0, 1], site_in, 100, cylinder.normal, site_type="in")
        mols_ex = pore.attach(mol, 0, [0, 1], site_ex, 20, normal, site_type="ex")

        ## Filling
        mols_in_fill = pore.fill_sites(site_in, cylinder.normal, site_type="in")
        mols_ex_fill = pore.fill_sites(site_ex, normal, site_type="ex")

        ## Storage
        pms.Store(pms.Molecule(name="pore_cylinder_siloxane", inp=mols_siloxane), "output").gro()
        pms.Store(pms.Molecule(name="pore_cylinder_in", inp=mols_in), "output").gro()
        pms.Store(pms.Molecule(name="pore_cylinder_ex", inp=mols_ex), "output").gro()
        pms.Store(pms.Molecule(name="pore_cylinder_in_fill", inp=mols_in_fill), "output").gro()
        pms.Store(pms.Molecule(name="pore_cylinder_ex_fill", inp=mols_ex_fill), "output").gro()

        # Delete atoms
        block.delete(matrix.bound(0))
        pms.Store(block, "output").gro()

        # Set reservoir
        pore.reservoir(5)
        self.assertEqual([round(x) for x in pore.get_box()], [6, 6, 17])

        # Output
        pore.set_name("pore_cylinder_full")
        pms.Store(pore, "output").gro(use_atom_names=True)

        pore.set_name("pore_cylinder_full_sort")
        sort_list = ["OM", "SI", "SLX", "SL", "SLG", "TMS", "TMSG"]
        pms.Store(pore, "output", sort_list=sort_list).gro(use_atom_names=True)
        pms.Store(pore, "output", sort_list=sort_list).pdb(use_atom_names=True)

        # Store test
        print()
        pms.Store(pore, "output", sort_list=sort_list[:-1])
        pms.Store(pore, "output", sort_list=sort_list).top()

        # Error test
        self.assertIsNone(pore.attach(mol, 0, [0, 1], site_in, 0, cylinder.normal, site_type="DOTA"))
        self.assertIsNone(pore.siloxane(site_in, 0, cylinder.normal, site_type="DOTA"))

        # Getter and Setter
        self.assertEqual(pore.get_block().get_name(), "pore_cylinder_block")
        self.assertEqual(len(pore.get_site_dict()), 3)
        self.assertEqual(pore.get_num_in_ex(), 23)

    def test_pore_exterior(self):
        # self.skipTest("Temporary")

        # x-axis
        pattern = pms.BetaCristobalit()
        block = pattern.generate([2, 2, 2], "x")
        block.set_name("pattern_beta_cbt_ex_x")
        dice = pms.Dice(block, 0.2, True)
        bonds = dice.find_bond(None, ["Si", "O"], [0.155-1e-2, 0.155+1e-2])
        matrix = pms.Matrix(bonds)
        pore = pms.Pore(block, matrix)
        pore.prepare()
        pore.exterior(pattern.get_gap())
        pore.sites()
        pms.Store(block, "output").gro()

        # y-axis
        pattern = pms.BetaCristobalit()
        block = pattern.generate([2, 2, 2], "y")
        block.set_name("pattern_beta_cbt_ex_y")
        dice = pms.Dice(block, 0.2, True)
        bonds = dice.find_bond(None, ["Si", "O"], [0.155-1e-2, 0.155+1e-2])
        matrix = pms.Matrix(bonds)
        pore = pms.Pore(block, matrix)
        pore.prepare()
        pore.exterior(pattern.get_gap())
        pore.sites()
        pms.Store(block, "output").gro()

        # z-axis
        pattern = pms.BetaCristobalit()
        block = pattern.generate([2, 2, 2], "z")
        block.set_name("pattern_beta_cbt_ex_z")
        dice = pms.Dice(block, 0.2, True)
        bonds = dice.find_bond(None, ["Si", "O"], [0.155-1e-2, 0.155+1e-2])
        matrix = pms.Matrix(bonds)
        pore = pms.Pore(block, matrix)
        pore.prepare()
        pore.exterior(pattern.get_gap())
        pore.sites()
        pms.Store(block, "output").gro()

        # Amorph
        pattern = pms.BetaCristobalit()
        pattern.generate([2, 2, 2], "z")

        pattern._structure = pms.Molecule(inp="data/amorph.gro")
        pattern._size = [2.014, 1.751, 2.468]

        block = pattern.get_block()
        block.set_name("pattern_beta_cbt_ex_amoprh")

        dice = pms.Dice(block, 0.4, True)
        matrix = pms.Matrix(dice.find_parallel(None, ["Si", "O"], [0.160-0.02, 0.160+0.02]))

        connect = matrix.get_matrix()
        matrix.split(57790, 2524)

        pore = pms.Pore(block, matrix)
        pore.prepare()
        pore.exterior(pattern.get_gap())
        pore.sites()

        pms.Store(block, "output").gro()

    def test_pore_cylinder(self):
        # self.skipTest("Temporary")

        # Empty pore
        pore = pms.PoreCylinder([4, 4, 4], 2, 0)
        pore.finalize()

        # Filled pore
        pore = pms.PoreCylinder([6, 6, 6], 4, 5, [5, 5])

        ## Attachement
        pore.attach_special(pms.gen.tms(),  0, [0, 1], 5)
        pore.attach_special(pms.gen.tms(),  0, [0, 1], 3, symmetry="mirror")

        tms2 = pms.gen.tms()
        tms2.set_short("TMS2")

        pore.attach(tms2, 0, [0, 1], 10, "in", trials=10, inp="percent")
        pore.attach(tms2, 0, [0, 1], 1, "in", trials=10, inp="molar")
        pore.attach(tms2, 0, [0, 1], 0.1, "ex", trials=10, inp="molar")

        # Special cases
        print()
        self.assertIsNone(pore.attach(pms.gen.tms(), 0, [0, 1], 100, site_type="DOTA"))
        self.assertIsNone(pore.attach(pms.gen.tms(), 0, [0, 1], 100, "in", inp="DOTA"))
        self.assertIsNone(pore.attach(pms.gen.tms(), 0, [0, 1], 100, pos_list=[[1, 3, 3], [7, 4, 2]]))
        self.assertIsNone(pore.attach_special(pms.gen.tms(),  0, [0, 1], 3, symmetry="DOTA"))

        # Finalize
        pore.finalize()
        pore.store("output/cylinder/")
        print(pore.table())

        ## Properties
        self.assertEqual(round(pore.diameter()), 4)
        self.assertEqual([round(x, 4) for x in pore.centroid()], [3.0773, 3.0934, 3.255])
        self.assertEqual(round(pore.roughness()["in"], 1), 0.1)
        self.assertEqual(round(pore.roughness()["ex"], 1), 0.0)
        self.assertEqual(round(pore.volume()), 80)
        self.assertEqual({key: round(item) for key, item in pore.surface().items()}, {'in': 79, 'ex': 49})
        self.assertEqual(pore.shape(), "CYLINDER")

    def test_pore_slit(self):
        # self.skipTest("Temporary")

        # Empty pore
        pore = pms.PoreSlit([4, 4, 4], 2)
        pore.finalize()

        # Filled pore
        pore = pms.PoreSlit([6, 6, 6], 3, 5, [5, 5])

        ## Attachement
        pore.attach_special(pms.gen.tms(),  0, [0, 1], 5)
        pore.attach_special(pms.gen.tms(),  0, [0, 1], 3, symmetry="mirror")

        tms2 = pms.gen.tms()
        tms2.set_short("TMS2")

        pore.attach(tms2, 0, [0, 1], 10, "in", trials=10, inp="percent")
        pore.attach(tms2, 0, [0, 1], 1, "in", trials=10, inp="molar")
        pore.attach(tms2, 0, [0, 1], 0.1, "ex", trials=10, inp="molar")

        # Special cases
        print()
        self.assertIsNone(pore.attach(pms.gen.tms(), 0, [0, 1], 100, site_type="DOTA"))
        self.assertIsNone(pore.attach(pms.gen.tms(), 0, [0, 1], 100, "in", inp="DOTA"))
        self.assertIsNone(pore.attach_special(pms.gen.tms(),  0, [0, 1], 3, symmetry="DOTA"))

        # Finalize
        pore.finalize()
        pore.store("output/slit/")
        print(pore.table())

        ## Properties
        self.assertEqual(round(pore.height()), 3)
        self.assertEqual([round(x, 4) for x in pore.centroid()], [3.0773, 3.0934, 3.255])
        self.assertEqual(round(pore.roughness()["in"], 1), 0.1)
        self.assertEqual(round(pore.roughness()["ex"], 1), 0.0)
        self.assertEqual(round(pore.volume()), 114)
        self.assertEqual(round(pore.surface()["in"]), 75)
        self.assertEqual(pore.shape(), "SLIT")

    def test_pore_capsule(self):
        # self.skipTest("Temporary")

        # Empty pore
        pore = pms.PoreCapsule([4, 4, 4], 1, 1, 0)
        pore.finalize()

        # Filled pore
        pore = pms.PoreCapsule([6, 6, 10], 4, 2, 5, [5, 5])

        ## Attachement
        # pore.attach_special(pms.gen.tms(),  0, [0, 1], 5)
        # pore.attach_special(pms.gen.tms(),  0, [0, 1], 3, symmetry="mirror")

        tms2 = pms.gen.tms()
        tms2.set_short("TMS2")

        pore.attach(tms2, 0, [0, 1], 10, "in", trials=10, inp="percent")
        pore.attach(tms2, 0, [0, 1], 1, "in", trials=10, inp="molar")
        pore.attach(tms2, 0, [0, 1], 0.1, "ex", trials=10, inp="molar")

        # Special cases
        print()
        self.assertIsNone(pore.attach(pms.gen.tms(), 0, [0, 1], 100, site_type="DOTA"))
        self.assertIsNone(pore.attach(pms.gen.tms(), 0, [0, 1], 100, "in", inp="DOTA"))
        # elf.assertIsNone(pore.attach_special(pms.gen.tms(),  0, [0, 1], 3, symmetry="DOTA"))

        # Finalize
        pore.finalize()
        pore.store("output/capsule/")
        print(pore.table())

        # Properties
        self.assertEqual(round(pore.diameter()), 4)
        self.assertEqual([round(x, 4) for x in pore.centroid()["block"]], [3.0775, 3.0935, 5.115])
        self.assertEqual(round(pore.roughness()["in"], 1), 0.1)
        self.assertEqual(round(pore.roughness()["ex"], 1), 0.0)
        self.assertEqual(round(pore.volume()), 86)
        self.assertEqual({key: round(item) for key, item in pore.surface().items()}, {'in': 102, 'ex': 49})
        self.assertEqual(pore.shape(), "CAPSULE")

    def test_pore_cylinder_amorph(self):
        # self.skipTest("Temporary")

        # Empty pore
        pore = pms.PoreAmorphCylinder(2, 0)
        pore.finalize()

        # Filled pore
        pore = pms.PoreAmorphCylinder(4, 5, [2, 2])

        ## Attachement
        pore.attach_special(pms.gen.tms(),  0, [0, 1], 5)
        pore.attach_special(pms.gen.tms(),  0, [0, 1], 3, symmetry="mirror")

        tms2 = pms.gen.tms()
        tms2.set_short("TMS2")

        pore.attach(tms2, 0, [0, 1], 10, "in", trials=10, inp="percent")
        pore.attach(tms2, 0, [0, 1], 1, "in", trials=10, inp="molar")
        pore.attach(tms2, 0, [0, 1], 0.1, "ex", trials=10, inp="molar")

        # Special cases
        print()
        self.assertIsNone(pore.attach(pms.gen.tms(), 0, [0, 1], 100, site_type="DOTA"))
        self.assertIsNone(pore.attach(pms.gen.tms(), 0, [0, 1], 100, "in", inp="DOTA"))
        self.assertIsNone(pore.attach(pms.gen.tms(), 0, [0, 1], 100, pos_list=[[1, 3, 3], [7, 4, 2]]))
        self.assertIsNone(pore.attach_special(pms.gen.tms(),  0, [0, 1], 3, symmetry="DOTA"))

        # Finalize
        pore.finalize()
        pore.store("output/cylinder_amorph/")
        print(pore.table())

        ## Properties
        self.assertEqual(round(pore.diameter()), 4)
        self.assertEqual([round(x, 4) for x in pore.centroid()], [4.923, 4.9651, 5.118])
        self.assertEqual(round(pore.roughness()["in"], 1), 0.1)
        self.assertEqual(round(pore.roughness()["ex"], 1), 0.2)
        self.assertEqual(round(pore.volume()), 122)
        self.assertEqual({key: round(item) for key, item in pore.surface().items()}, {'in': 121, 'ex': 159})
        self.assertEqual(pore.shape(), "CYLINDER")


if __name__ == '__main__':
    unittest.main(verbosity=2)
