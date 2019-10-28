import os
import sys

sys.path.append(os.path.abspath('../porems'))

from molecule import Molecule

# Generate molecules
molGro = Molecule(inp="data/benzene.gro")
molPdb = Molecule(inp="data/benzene.pdb")
molMol = Molecule(inp="data/benzene.mol2")

# _read
for i in range(4):
    print(molGro.getData()[i])
    print(molPdb.getData()[i])
    print(molMol.getData()[i])

# _concat
for x in Molecule([molGro,molGro]).getData(): print(x)
