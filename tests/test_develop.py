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


    ###########
    # Pattern #
    ###########
    def test_pattern_beta_cristobalit(self):
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

        # Exterior
        beta_cristobalit = pms.BetaCristobalit()
        beta_cristobalit.generate([2, 2, 2], "x")
        beta_cristobalit.get_block().set_name("pattern_beta_cbt_ex_x")
        beta_cristobalit.exterior()
        self.assertEqual(beta_cristobalit.get_block().get_num(), 600)
        pms.Store(beta_cristobalit.get_block(), "output").gro()

        beta_cristobalit = pms.BetaCristobalit()
        beta_cristobalit.generate([2, 2, 2], "y")
        beta_cristobalit.get_block().set_name("pattern_beta_cbt_ex_y")
        beta_cristobalit.exterior()
        self.assertEqual(beta_cristobalit.get_block().get_num(), 608)
        pms.Store(beta_cristobalit.get_block(), "output").gro()

        beta_cristobalit = pms.BetaCristobalit()
        beta_cristobalit.generate([2, 2, 2], "z")
        beta_cristobalit.get_block().set_name("pattern_beta_cbt_ex_z")
        beta_cristobalit.exterior()
        self.assertEqual(beta_cristobalit.get_block().get_num(), 592)
        pms.Store(beta_cristobalit.get_block(), "output").gro()

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


    #########
    # Shape #
    #########
    def test_shape_hyperboloid(self):
        # self.skipTest("Temporary")

        block = pms.BetaCristobalit().generate([6, 6, 6], "z")
        block.set_name("shape_hyperboloid")
        dice = pms.Dice(block, 0.4, True)
        matrix = pms.Matrix(dice.find_parallel(None, ["Si", "O"], 0.155, 1e-2))
        centroid = block.centroid()
        central = pms.geom.unit(pms.geom.rotate([0, 0, 1], [1, 0, 0], 45, True))
        central = [0, 0, 1]

        hyperboloid = pms.Hyperboloid({"centroid": centroid, "central": central, "length": 2, "base": 2, "skirt": 2})

        # Properties
        # self.assertEqual(round(cylinder.volume(), 4), 37.6991)
        # self.assertEqual(round(cylinder.surface(), 4), 37.6991)

        # Test vector
        vec = [3.6086, 4.4076, 0.2065]

        # Surface
        # self.assertEqual([round(x[0][20], 4) for x in cylinder.surf(num=100)], vec)
        # self.assertEqual([round(x[0][20], 4) for x in cylinder.rim(0, num=100)], vec)

        # Normal
        # self.assertEqual([round(x, 4) for x in cylinder.convert([0, 0, 0], False)], [3.0147, 3.0572, 1.5569])
        # self.assertEqual([round(x, 4) for x in cylinder.normal(vec)], [0.5939, 2.9704, 0.0000])

        # Positioning
        del_list = [atom_id for atom_id, atom in enumerate(block.get_atom_list()) if hyperboloid.is_in(atom.get_pos())]
        matrix.strip(del_list)
        block.delete(matrix.bound(0))
        # self.assertEqual(block.get_num(), 12650)

        # Store molecule
        pms.Store(block, "output").gro()

        # Plot surface
        plt.figure()
        hyperboloid.plot()#vec=[3.17290646, 4.50630614, 0.22183271])
        plt.show()


    ########
    # Pore #
    ########
    def test_pore_cylinder(self):
        self.skipTest("Temporary")

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

        ## Properties
        self.assertEqual(round(pore.diameter()), 4)
        self.assertEqual([round(x, 4) for x in pore.centroid()], [3.0773, 3.0934, 3.255])
        self.assertEqual(round(pore.roughness(), 1), 0.1)
        # self.assertEqual(round(pore.volume()), 81)
        # self.assertEqual({key: round(item) for key, item in pore.surface().items()}, {'in': 81, 'ex': 49})

        print(pore.table())


if __name__ == '__main__':
    unittest.main(verbosity=2)
