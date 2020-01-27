################################################################################
# Store Class                                                                  #
#                                                                              #
"""Contains function for creating simulation files."""
################################################################################


import os
import shutil

from decimal import Decimal

import porems.utils as utils

from porems.molecule import Molecule
from porems.pore import Pore


class Store:
    """This class converts a molecule object to a structure file of different
    possible formats. Currently only **PDB** and **GRO** are fully supported.

    For pore objects two additional methods are available for generating the
    main topology file and a topology file containing grid molecule
    parameters and charges. Currently only the GROMACS format is supported.

    Furthermore, there is an automized routine for generating topologies
    with Antechamber where job-file for this tool are created.

    Parameters
    ----------
    molecule : Molecule
        Molecule object
    link : string, optional
        Folder link for output

    Examples
    --------
    Assuming a pore or molecule object have already been created, the structure
    files can be generated like following examples

    .. code-block:: python

        from porems.store import Store

        Store(pore).gro()
        Store(mol).pdb()
    """
    def __init__(self, molecule, link="./"):
        # Get molecule properties
        self._dim = 3
        self._link = link if link[-1] == "/" else link+"/"
        self._mol = molecule
        self._mols = molecule.get_mol_list()
        self._name = molecule.get_name()
        self._short = molecule.get_short()
        self._box = Molecule(inp=self._mols).get_box() if molecule.get_box_c() == None else molecule.get_box_c()

        # Create output folder
        utils.mkdirp(link)


    ################################
    # Public Methods - Antechamber #
    ################################
    def job(self, is_master=False):
        """Create job file to run with Antechamber. A shell file and a tleap
        file are created with all necessary commands to create a topology from a
        pdb file. For the conversion to GROMACS file format the python package
        **ParmED** needs to be installed.

        Parameters
        ----------
        is_master : bool, optional
            True if the jobfile call should be added to a master run file
            (practical for running multiple topology generations)
        """
        # Initialize
        name = self._name.lower()
        short = self._short
        link = self._link

        # Template directory
        package_dir = os.path.split(__file__)[0]+"/"

        # Job file
        file_in = package_dir+"templates/antechamber.job"
        file_out = link+name+".job"
        shutil.copy(file_in, file_out)

        utils.replace(file_out, "MOLNAME", name)

        # Tleap file
        file_in = package_dir+"templates/antechamber.tleap"
        file_out = link+name+".tleap"
        shutil.copy(file_in, file_out)

        utils.replace(file_out, "MOLSHORTLOWER", name.lower())
        utils.replace(file_out, "MOLSHORT", short)
        utils.replace(file_out, "MOLNAME", name)

        # Add to master run
        if is_master:
            fileMaster = open(link+"run.job", "a")
            fileMaster.write("cd "+name+" #\n")
            fileMaster.write("sh "+name+".job #\n")
            fileMaster.write("cd .. #\n")
            fileMaster.write("echo \"Finished "+name+"...\"\n")
            fileMaster.write("#\n")
            fileMaster.close()


    ##############################
    # Public Methods - Structure #
    ##############################
    def obj(self, name=None):
        """Save the molecule object using pickle.

        Parameters
        ----------
        name : None, string, optional
            Filename
        """
        # Initialize
        link = self._link
        link += self._name+".obj" if name is None else name

        # Save object
        utils.save(self._mol, link)

    def pdb(self, name=None):
        """Generate the structure file for the defined molecule in the **PDB**
        format.

        Parameters
        ----------
        name : None, string, optional
            Filename
        """
        # Initialize
        dim = self._dim
        mols = self._mols
        link = self._link
        link += self._name+".pdb" if name is None else name

        # Open file
        with open(link, "w") as file_out:
            # Atoms
            num_a = 1
            num_m = 1
            for mol_id in range(len(mols)):
                data = mols[mol_id].get_data()
                short = mols[mol_id].get_short()
                ids = {x: 1 for x in data[dim]}
                for i in range(len(data[0])):
                    atom = data[dim][i]
                    out_string = "HETATM"                        #  1- 6 (6)    Record name
                    out_string += "%5i" % num_a                  #  7-11 (5)    Atom serial number
                    out_string += " "                            # 12    (1)    -
                    out_string += "%4s" % (atom+str(ids[atom]))  # 13-16 (4)    Atom name
                    # 17    (1)    Alternate location indicator
                    out_string += " "
                    out_string += "%3s" % short                  # 18-20 (3)    Residue name
                    out_string += " "                            # 21    (1)    -
                    out_string += "%1s" % "A"                    # 22    (1)    Chain identifier
                    # 23-26 (4)    Residue sequence number
                    out_string += "%4i" % (num_m)
                    # 27    (1)    Code for insertion of residues
                    out_string += " "
                    out_string += "   "                          # 28-30 (3)    -

                    for j in range(dim):                         # 31-54 (3*8)  Coordinates
                        out_string += "%8.3f" % (data[j][i]*10)

                    out_string += "%6.2f" % 1                    # 55-60 (6)    Occupancy
                    out_string += "%6.2f" % 0                    # 61-66 (6)    Temperature factor
                    out_string += "           "                  # 67-77 (11)   -
                    out_string += "%2s" % atom                   # 78-79 (2)    Element symbol
                    out_string += "  "                           # 80-81 (2)    Charge on the atom

                    file_out.write(out_string+"\n")

                    ids[atom] += 1
                    num_a += 1

                num_m = num_m+1 if num_m < 9999 else 1

            # End statement
            file_out.write("TER\nEND\n")

    def gro(self, name=None):
        """Generate the structure file for the defined molecule in the **GRO**
        format.

        Parameters
        ----------
        name : None, string, optional
            Filename
        """
        # Initialize
        dim = self._dim
        box = self._box
        mols = self._mols
        link = self._link
        link += self._name+".gro" if name is None else name

        # Open file
        with open(link, "w") as file_out:
            # Set title
            file_out.write("Molecule generated by the PoreMS package\n")

            # Number of atoms
            file_out.write("%5i" % sum([len(x.get_data()[0]) for x in mols])+"\n")

            # Atoms
            num_a = 1
            num_m = 1
            for mol_id in range(len(mols)):
                data = mols[mol_id].get_data()
                short = mols[mol_id].get_short()
                ids = {x: 1 for x in data[dim]}
                for i in range(len(data[0])):
                    atom = data[dim][i]
                    atom_name = atom+str(ids[atom])

                    out_string = "%5i" % num_m       # 1- 5 (5)    Residue number
                    out_string += "%-5s" % short     # 6-10 (5)    Residue short name
                    out_string += "%5s" % atom_name  # 11-15 (5)    Atom name
                    out_string += "%5i" % num_a      # 16-20 (5)    Atom number
                    for j in range(dim):             # 21-44 (3*8)  Coordinates
                        out_string += "%8.3f" % data[j][i]

                    file_out.write(out_string+"\n")

                    num_a = num_a+1 if num_a < 99999 else 0
                    ids[atom] += 1

                num_m = num_m+1 if num_m < 99999 else 0

            # Box
            out_string = ""
            for i in range(dim):
                out_string += "%.3f" % box[i]
                out_string += " " if i < dim-1 else "\n"

            file_out.write(out_string)

    def xyz(self, name=None):
        """Generate the structure file for the defined molecule in the **XYZ**
        format for running qm-simulations.

        Parameters
        ----------
        name : string, None, optional
            Filename
        """
        # Initialize
        mols = self._mols
        link = self._link
        link += self._name+".xyz" if name is None else name

        # Open output file and set title
        with open(link, "w") as file_out:
            file_out.write(str(sum([len(x.get_data()[0]) for x in mols]))+"\n"+"Energy = \n")

            # Atoms
            for mol_id in range(len(mols)):
                data = mols[mol_id].get_data()
                for i in range(len(data[0])):
                    out_string = "%-2s" % data[-1][i]  # 1- 2 (2)     Atom name

                    for j in range(self._dim):       # 3-41 (3*13)  Coordinates
                        out_string += "%13.7f" % (data[j][i]*10)

                    file_out.write(out_string+"\n")


    ############
    # Topology #
    ############
    def top(self, name=None):
        """Store the **topology** file for a pore. A top file is created
        containing the itp-include for all molecules and the count of the
        different groups of the pore.

        Parameters
        ----------
        name : None, string, optional
            Filename
        """
        # Initialize
        mols = self._mols
        link = self._link
        link += self._name+".top" if name is None else name

        # Get unique molecules
        unique_mols = []
        for mol in mols:
            if not mol.get_name() in unique_mols:
                unique_mols.append(mol.get_name())

        # Open file
        with open(link, "w") as file_out:
            # Store header
            file_out.write("[ defaults ]\n")
            file_out.write("; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ\n")
            file_out.write("1               2               yes             0.5     0.833333\n")
            file_out.write("\n")
            file_out.write("#include \"grid.itp\"\n")

            for mol_name in unique_mols:
                if mol_name not in ["si", "om", "ox", "sl", "slg", "slx"]:
                    file_out.write("#include \""+mol_name+".itp\"\n")

            file_out.write("\n")
            file_out.write("[ system ]\n")
            file_out.write("A pore used in molecular simulation plus reservoir\n")
            file_out.write("\n")
            file_out.write("[ molecules ]\n")

            # Atoms
            counter = 1
            for i in range(1, len(mols)):
                if mols[i].get_name() == mols[i-1].get_name():
                    counter += 1
                else:
                    file_out.write(mols[i-1].get_short()+" "+str(counter)+"\n")
                    counter = 1

            file_out.write(mols[-1].get_short()+" "+str(counter)+"\n")

    def grid(self, name=None):
        """Store the **grid.itp** file containing the necessary parameters and
        charges of the grid molecules.

        Parameters
        ----------
        name : None, string, optional
            Filename
        """
        # Initialize
        link = self._link
        link += "grid.itp" if name is None else name

        # Calculate excess charge
        charges = self._mol.get_q()

        # Copy grid file
        utils.copy(os.path.split(__file__)[0]+"/templates/grid.itp", link)

        # Replace charges
        utils.replace(link, "CHARGEOX", "%8.6f" % charges["ox"])
        utils.replace(link, "CHARGEO", "%8.6f" % charges["om"])
        utils.replace(link, "CHARGESI", "%8.6f" % charges["si"])
