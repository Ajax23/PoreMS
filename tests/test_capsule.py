import os
import sys

import porems as pms


########
# Cube #
########
def shape_capsule():
    pattern = pms.BetaCristobalit()
    block = pattern.generate([6, 6, 12], "z")
    block.set_name("shape_capsule")
    dice = pms.Dice(block, 0.4, True)
    matrix = pms.Matrix(dice.find_parallel(None, ["Si", "O"], 0.155, 10e-2))
    central = pms.geom.unit(pms.geom.rotate([0, 0, 1], [1, 0, 0], 0, True))
    centroid = block.centroid()
    centroid_cyl_l = centroid[:2]+[0]
    centroid_cyl_r = centroid[:2]+[pattern.get_size()[2]]
    centroid_sph_l = centroid[:2]+[3]
    centroid_sph_r = centroid[:2]+[pattern.get_size()[2]-3]

    cylinder_l = pms.Cylinder({"centroid": centroid_cyl_l, "central": central, "length": 6, "diameter": 4})
    cylinder_r = pms.Cylinder({"centroid": centroid_cyl_r, "central": central, "length": 6, "diameter": 4})
    sphere_l = pms.Sphere({"centroid": centroid_sph_l, "central": central, "diameter": 4})
    sphere_r = pms.Sphere({"centroid": centroid_sph_r, "central": central, "diameter": 4})

    del_list = []
    del_list.extend([atom_id for atom_id, atom in enumerate(block.get_atom_list()) if cylinder_l.is_in(atom.get_pos())])
    del_list.extend([atom_id for atom_id, atom in enumerate(block.get_atom_list()) if cylinder_r.is_in(atom.get_pos())])
    del_list.extend([atom_id for atom_id, atom in enumerate(block.get_atom_list()) if sphere_l.is_in(atom.get_pos())])
    del_list.extend([atom_id for atom_id, atom in enumerate(block.get_atom_list()) if sphere_r.is_in(atom.get_pos())])
    matrix.strip(del_list)

    block.delete(matrix.bound(0))
    pms.Store(block, "output").gro()

def pore_capsule():
    orient = "z"
    pattern = pms.BetaCristobalit()
    pattern.generate([6, 6, 12], orient)
    pattern.exterior()

    block = pattern.get_block()
    block.set_name("pore_capsule")

    dice = pms.Dice(block, 0.4, True)
    matrix = pms.Matrix(dice.find_parallel(None, ["Si", "O"], 0.155, 10e-2))
    oxygen_out = matrix.bound(1)

    centroid = block.centroid()
    central = pms.geom.unit(pms.geom.rotate([0, 0, 1], [1, 0, 0], 0, True))
    centroid_cyl_l = centroid[:2]+[0]
    centroid_cyl_r = centroid[:2]+[pattern.get_size()[2]]
    centroid_sph_l = centroid[:2]+[3]
    centroid_sph_r = centroid[:2]+[pattern.get_size()[2]-3]

    cylinder_l = pms.Cylinder({"centroid": centroid_cyl_l, "central": central, "length": 6, "diameter": 4})
    cylinder_r = pms.Cylinder({"centroid": centroid_cyl_r, "central": central, "length": 6, "diameter": 4})
    sphere_l = pms.Sphere({"centroid": centroid_sph_l, "central": central, "diameter": 4})
    sphere_r = pms.Sphere({"centroid": centroid_sph_r, "central": central, "diameter": 4})

    def normal_in(pos):
        if pos[2] <= 3:
            return cylinder_l.normal(pos)
        elif pos[2] > 3 and pos[2] < centroid[2]:
            return [x if i<2 else -x for i, x in enumerate(sphere_l.normal(pos))]
        elif pos[2] > centroid[2] and pos[2] < pattern.get_size()[2]-3:
            return [x if i<2 else -x for i, x in enumerate(sphere_r.normal(pos))]
        elif pos[2] >= pattern.get_size()[2]-3:
            return cylinder_r.normal(pos)

    def normal_ex(pos):
        return [0, 0, -1] if pos[2] < centroid[2] else [0, 0, 1]

    del_list = []
    del_list.extend([atom_id for atom_id, atom in enumerate(block.get_atom_list()) if cylinder_l.is_in(atom.get_pos())])
    del_list.extend([atom_id for atom_id, atom in enumerate(block.get_atom_list()) if cylinder_r.is_in(atom.get_pos())])
    del_list.extend([atom_id for atom_id, atom in enumerate(block.get_atom_list()) if sphere_l.is_in(atom.get_pos())])
    del_list.extend([atom_id for atom_id, atom in enumerate(block.get_atom_list()) if sphere_r.is_in(atom.get_pos())])
    matrix.strip(del_list)

    # Process surface
    pore = pms.Pore(block, matrix)
    pore.prepare()
    pore.sites(oxygen_out)
    site_list = pore.get_sites()
    site_in = [site_key for site_key, site_val in site_list.items() if site_val["type"]=="in"]
    site_ex = [site_key for site_key, site_val in site_list.items() if site_val["type"]=="ex"]

    # Attachment
    mol = pms.gen.tms()

    mols_in = pore.attach(mol, 0, [0, 1], site_in, 100, normal_in)
    mols_ex = pore.attach(mol, 0, [0, 1], site_ex, 20, normal_ex)
    mols_in_fill = pore.fill_sites(site_in, normal_in)
    mols_ex_fill = pore.fill_sites(site_ex, normal_ex)

    pms.Store(pms.Molecule(name="pore_capsule_in", inp=mols_in), "output").gro()
    pms.Store(pms.Molecule(name="pore_capsule_ex", inp=mols_ex), "output").gro()
    pms.Store(pms.Molecule(name="pore_capsule_in_fill", inp=mols_in_fill), "output").gro()
    pms.Store(pms.Molecule(name="pore_capsule_ex_fill", inp=mols_ex_fill), "output").gro()

    # Objectify
    mol_obj = pore.objectify()
    pms.Store(pms.Molecule(name="pore_capsule_grid", inp=mol_obj), "output").gro(use_atom_names=True)

    # Delete atoms
    block.delete(matrix.bound(0))
    pms.Store(block, "output").gro()

    # Output
    pore.set_box(block.get_box())

    pore.set_name("pore_capsule_full")
    pms.Store(pore, "output").gro(use_atom_names=True)

    pore.set_name("pore_capsule_full_sort")
    sort_list = ["OM", "SI", "SL", "SLG", "TMS", "TMSG"]
    pms.Store(pore, "output", sort_list=sort_list).gro(use_atom_names=True)


if __name__ == '__main__':
    shape_capsule()
    pore_capsule()
