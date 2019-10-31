import os
import sys

# Install package
os.system("pip install ../. &> /dev/null")
print("Finished inistalling package...")

# import package
from porems import *

# Generate molecules
mol_gro = Molecule(inp="data/benzene.gro")
mol_pdb = Molecule(inp="data/benzene.pdb")
mol_mol = Molecule(inp="data/benzene.mol2")


############
# Molecule #
############
# # _read
# for i in range(4):
#     print(mol_gro.get_data()[i])
#     print(mol_pdb.get_data()[i])
#     print(mol_mol.get_data()[i])
#
# # _concat
# for x in Molecule([mol_gro,mol_gro]).get_data(): print(x)

# set_masses
# print(mol_gro.get_masses())

#########
# Write #
#########
# Write(mol_gro,"output").job()
# Write(mol_gro,"output").gro()
# Write(mol_gro,"output").pdb()
# Write(mol_gro,"output").xyz()


########
# Pore #
########
pore = Pore(size=[10, 10, 10], diam=6, drill="z", res=5.5, is_time=True)
pore.set_name("pore")

# pore.couple([C18(),C18()],    [[0,7,8],[0,7,8]],[30.4,39])
# pore.couple([Silyl(),Silyl()],[[0,1,2],[0,1,2]],[54-30.4,54-39])
#
# pore.setGrid([C18(),Silyl()])

pore.finish(is_mol=True, is_props=False)

Write(pore, "output").gro()
Write(pore, "output").top()
Write(pore, "output").grid()
