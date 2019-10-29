import os
import sys

# Install package
os.system("pip install ../. &> /dev/null")
print("Finished inistalling package...")

# import package
from porems import *


# Generate molecules
molGro = Molecule(inp="data/benzene.gro")
molPdb = Molecule(inp="data/benzene.pdb")
molMol = Molecule(inp="data/benzene.mol2")


############
# Molecule #
############
# # _read
# for i in range(4):
#     print(molGro.get_data()[i])
#     print(molPdb.get_data()[i])
#     print(molMol.get_data()[i])
#
# # _concat
# for x in Molecule([molGro,molGro]).get_data(): print(x)


#########
# Write #
#########
# Write(molGro,"output").job()
# Write(molGro,"output").gro()
# Write(molGro,"output").pdb()
# Write(molGro,"output").xyz()


########
# Pore #
########
pore = Pore(size=[10,10,10],diam=6,drill="z",res=5.5,is_time=True)
pore.set_name("pore")

# pore.couple([C18(),C18()],    [[0,7,8],[0,7,8]],[30.4,39])
# pore.couple([Silyl(),Silyl()],[[0,1,2],[0,1,2]],[54-30.4,54-39])
#
# pore.setGrid([C18(),Silyl()])

pore.finish(is_mol=True,is_props=False)

Write(pore,"output").gro()
Write(pore,"output").top()
Write(pore,"output").grid()
