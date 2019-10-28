################################################################################
# Write Class                                                                  #
#                                                                              #
"""Contains function for creating simululation files."""
################################################################################


import os
import shutil

from decimal import Decimal

import porems.utils as utils

from porems.molecule import Molecule
from porems.pore     import Pore


class Write:
    """This class converts a molecule object to a text file of different
    possible formats. Currently only **PDB** and **GRO** are fully supported.

    Furthermore there is an automized rountine for generating topologies
    with Antechamber where job-file for this tool are created.

    Parameters
    ----------
    molecule : Molecule
        Molecule object
    link : str
        Folder link for output
    """
    def __init__(self,molecule,link="./"):
        # Get molecule properties
        self._dim   = 3
        self._link  = link if link[-1]=="/" else link+"/"
        self._mols  = molecule.get_write()
        self._name  = molecule.get_name()
        self._short = molecule.get_short()
        self._box   = Molecule(self._mols).get_box() if molecule.get_box_c()==None else molecule.get_box_c()

        # Create output folder
        utils.mkdirp(link)


    ################################
    # Public Methods - Antechamber #
    ################################
    def job(self,is_master=False):
        """Create job file to run with Antechamber. A shell file and a tleap file
        are created with all necesarry commands to create a topology from a pdb file.
        For the conversion to gromacs file format the python package **ParmED**
        needs to be installed.

        Parameters
        ----------
        is_master : bool
            True if the jobfile call should be added to a master run file
            (pracitcall for running multiple topology genrations)
        """
        # Initialize
        name  = self._name.lower()
        short = self._short
        link  = self._link

        # Template directory
        package_dir = os.path.split(__file__)[0]+"/"

        # Job file
        file_in  = package_dir+"templates/antechamber.job"
        file_out = link+name+".job"
        shutil.copy(file_in,file_out)

        utils.replace(file_out,"MOLNAME",name)

        # Tleap file
        file_in  = package_dir+"templates/antechamber.tleap"
        file_out = link+name+".tleap"
        shutil.copy(file_in,file_out)

        utils.replace(file_out,"MOLSHORTLOWER",name.lower())
        utils.replace(file_out,"MOLSHORT",short)
        utils.replace(file_out,"MOLNAME",name)

        # Add to master run
        if is_master:
            fileMaster = open(link+"run.job","a")
            fileMaster.write("cd "+name+" #\n")
            fileMaster.write("sh "+name+".job #\n")
            fileMaster.write("cd .. #\n")
            fileMaster.write("echo \"Finished "+name+"...\"\n")
            fileMaster.write("#\n")
            fileMaster.close()


    ##############################
    # Public Methods - Structure #
    ##############################
    def pdb(self,name=None,link=None):
        """Generate the structure file for the defined molecule in the **PDB** format.
        If parameter *link* is given, parameter *name* will be ignored.

        Parameters
        ----------
        name : str, None
            Filename
        link : str, None
            Full link with filename
        """
        # Initialize
        dim   = self._dim
        mols  = self._mols
        if link is None:
            link  = self._link
            link += self._name+".pdb" if name is None else name

        # Open file
        with open(link,"w") as file_out:
            # Atoms
            num_a = 1
            num_m = 1
            for mol_id in range(len(mols)):
                data  = mols[mol_id].get_data()
                short = mols[mol_id].get_short()
                ids   = {x: 1 for x in data[dim]}
                for i in range(len(data[0])):
                    atom        = data[dim][i]
                    out_string  = "HETATM"                     #  1- 6 (6)    Record name
                    out_string += "%5i"%num_a                  #  7-11 (5)    Atom serial number
                    out_string += " "                          # 12    (1)    -
                    out_string += "%4s"%(atom+str(ids[atom]))  # 13-16 (4)    Atom name
                    out_string += " "                          # 17    (1)    Alternate location indicator
                    out_string += "%3s"%short                  # 18-20 (3)    Residue name
                    out_string += " "                          # 21    (1)    -
                    out_string += "%1s"%"A"                    # 22    (1)    Chain identifier
                    out_string += "%4i"%(num_m)                # 23-26 (4)    Residue sequence number
                    out_string += " "                          # 27    (1)    Code for insertion of residues
                    out_string += "   "                        # 28-30 (3)    -

                    for j in range(dim):                       # 31-54 (3*8)  Corrdinates
                        out_string += "%8.3f"%(data[j][i]*10)

                    out_string += "%6.2f"%1                    # 55-60 (6)    Occupancy
                    out_string += "%6.2f"%0                    # 61-66 (6)    Temperatur factor
                    out_string += "           "                # 67-77 (11)   -
                    out_string += "%2s"%atom                   # 78-79 (2)    Element symbol
                    out_string += "  "                         # 80-81 (2)    Charge on the atom

                    file_out.write(out_string+"\n")

                    ids[atom] += 1
                    num_a     += 1

                num_m = num_m+1 if num_m<9999 else 1

            # End statement
            file_out.write("TER\nEND\n")

    def gro(self,name=None,link=None):
        """Generate the structure file for the defined molecule in the **GRO** format.
        If parameter *link* is given, parameter *name* will be ignored.

        Parameters
        ----------
        name : str, None
            Filename
        link : str, None
            Full link with filename
        """
        # Initialize
        dim  = self._dim
        box  = self._box
        mols = self._mols
        if link is None:
            link  = self._link
            link += self._name+".gro" if name is None else name

        # Open file
        with open(link,"w") as file_out:
            # Set title
            file_out.write("Molecule generated by the PoreMS package\n")

            # Number of atoms
            file_out.write("%5i"%sum([len(x.get_data()[0]) for x in mols])+"\n")

            # Atoms
            num_a = 1
            num_m = 1
            for mol_id in range(len(mols)):
                data  = mols[mol_id].get_data()
                short = mols[mol_id].get_short()
                ids   = {x: 1 for x in data[dim]}
                for i in range(len(data[0])):
                    atom        = data[dim][i]
                    atom_name   = atom+str(ids[atom])

                    out_string  = "%5i"% num_m      #  1- 5 (5)    Residue number
                    out_string += "%-5s"%short      #  6-10 (5)    Residue short name
                    out_string += "%5s"% atom_name  # 11-15 (5)    Atom name
                    out_string += "%5i"% num_a      # 16-20 (5)    Atom number
                    for j in range(dim):            # 21-44 (3*8)  Coordinates
                        out_string += "%8.3f"%data[j][i]

                    file_out.write(out_string+"\n")

                    num_a      = num_a+1 if num_a<99999 else 0
                    ids[atom] += 1

                num_m = num_m+1 if num_m<99999 else 0

            # Box
            out_string = ""
            for i in range(dim):
                out_string += "%.3f"%box[i]
                out_string += " " if i<dim-1 else "\n"

            file_out.write(out_string)

    def xyz(self,name=None,link=None):
        """Generate the structure file for the defined molecule in the xyz
        format for running qm-simulations.
        If parameter *link* is given, parameter *name* will be ignored.

        Parameters
        ----------
        name : str, None
            Filename
        link : str, None
            Full link with filename
        """
        # Initialize
        mols  = self._mols
        if link is None:
            link  = self._link
            link += self._name+".xyz" if name is None else name

        # Open output file and set title
        with open(link,"w") as file_out:
            file_out.write(str(sum([len(x.get_data()[0]) for x in mols]))+"\n"+"Energy = \n")

            # Atoms
            for mol_id in range(len(mols)):
                data = mols[mol_id].get_data()
                for i in range(len(data[0])):
                    out_string = "%-2s"%data[-1][i]  # 1- 2 (2)     Atom name

                    for j in range(self._dim):       # 3-41 (3*13)  Coordinates
                        out_string += "%13.7f"%(data[j][i]*10)

                    file_out.write(out_string+"\n")
