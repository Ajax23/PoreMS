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
    from porems import *
    from porems.essentials import Alkane
    from porems.essentials import Alcohol
    from porems.essentials import Ketone
    from porems.essentials import TMS

    # Generate molecules
    mol_gro = Molecule(inp="data/benzene.gro")
    mol_pdb = Molecule(inp="data/benzene.pdb")
    mol_mol = Molecule(inp="data/benzene.mol2")


    ############
    # Molecule #
    ############
    # _read
    for i in range(4):
        print(mol_gro.get_data()[i])
        print(mol_pdb.get_data()[i])
        print(mol_mol.get_data()[i])

    # _concat
    for x in Molecule(inp=[mol_gro, mol_gro]).get_data(): print(x)

    # set_masses
    print(mol_gro.get_masses())

    # table
    print(mol_gro.table())


    #########
    # Store #
    #########
    Store(mol_gro, "output").job()
    Store(mol_gro, "output").obj()
    Store(mol_gro, "output").gro()
    Store(mol_gro, "output").pdb()
    Store(mol_gro, "output").xyz()

    print(Molecule(inp="output/Molecule.obj").get_data())


    ##############
    # Essentials #
    ##############
    alkane = Alkane(10, "alkane")
    alcohol = Alcohol(10, "alcohol")
    ketone = Ketone(10,5, "ketone")
    tms = TMS(compress=30)

    Store(alkane, "output").gro()
    Store(alcohol, "output").gro()
    Store(ketone, "output").gro()
    Store(tms, "output").gro()


    ########
    # Pore #
    ########
    # Generation
    pore = Pore(size=[10, 10, 10], diam=6, drill="z", res=5, is_time=True)
    pore.set_name("pore")

    pore.siloxane(10, "num")
    pore.attach(TMS(compress=30), [0, 1], [1, 2], 0, 3, inp="percent")
    pore.attach(TMS(compress=30), [0, 1], [1, 2], 0, 0.67, inp="molar")
    pore.attach(TMS(compress=30), [0, 1], [1, 2], 1, 10, inp="num")

    pore.finalize()

    # Analysis
    props = pore.get_props()

    for prop in props:
        print(prop,":",props[prop])

    # Tables
    tables = pore.table_props()
    print(tables["props"])
    print(tables["alloc"])
    print(tables["time"])

    # Structure
    Store(pore, "output").obj()
    Store(pore, "output").gro()
    Store(pore, "output").top()
    Store(pore, "output").grid()
