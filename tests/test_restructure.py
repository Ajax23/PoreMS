import os
import sys

# Run test
if __name__ == "__main__":
    # Install package
    if sys.platform == "win32":
        os.system("pip install ../.")
    else:
        os.system("pip install ../. &> /dev/null")
    print("Finished inistalling package...")

    # Import package
    from porems.atom import Atom
    from porems.molecule2 import Molecule
    from porems.pore2 import Pore
    from porems.store2 import Store

    from porems.pattern import BetaCristobalit


    #########
    # Atoms #
    #########
    # print(Atom([0.0, 0.1, 0.2], "H", "LOL"))


    ############
    # Molecule #
    ############
    mol = Molecule(inp="data/benzene.gro")
    mol.set_name("molecule2")
    # mol_concat = Molecule(inp=[mol, mol])
    # mol_atom = Molecule(inp=mol.get_atom_list())
    # print(mol_atom)


    #########
    # Store #
    #########
    # Store(mol, "output").gro()
    # Store(mol, "output").pdb()
    # Store(mol, "output").xyz()


    ###########
    # Pattern #
    ###########
    # Store(BetaCristobalit().pattern(), "output").gro()

    ########
    # Pore #
    ########
    # pore = Pore([10, 10, 10], "z")
    # pore.generate()
    #
    # Store(pore.get_pore(), "output").gro()
