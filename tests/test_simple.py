import os
import sys

from porems.molecule import Molecule
from porems.write    import Write


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
