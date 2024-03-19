################################################################################
# Store Class                                                                  #
#                                                                              #
"""Contains function for creating simulation files."""
################################################################################


import os
import shutil

import porems.utils as utils
import porems.database as db

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
    inp : Molecule, Pore
        Molecule or Pore object
    link : string, optional
        Folder link for output
    sort_list : list, optional
        Optionally provide a sorting list in case a dictionary of molecules is
        given

    Examples
    --------
    Assuming a pore or molecule object have already been created, the structure
    files can be generated like following examples

    .. code-block:: python

        Store(mol).pdb()
        Store(pore, "output").gro("pore.gro")
    """
    def __init__(self, inp, link="./", sort_list=[]):
        # Initialize
        self._dim = 3
        self._link = link if link[-1] == "/" else link+"/"
        self._inp = inp

        # Process input
        if isinstance(inp, Molecule):
            self._mols = [self._inp]
        elif isinstance(inp, Pore):
            if sort_list:
                if sorted(sort_list) == sorted(list(inp.get_mol_dict().keys())):
                    self._mols = sum([inp.get_mol_dict()[x] for x in sort_list], [])
                    self._short_list = sort_list
                else:
                    print("Store: Sorting list does not contain all keys...")
                    return
            else:
                self._mols = sum([x for x in inp.get_mol_dict().values()], [])
                self._short_list = list(x for x in inp.get_mol_dict().keys())
        else:
            print("Store: Unsupported input type...")
            return

        # Get properties after input checking type
        self._name = inp.get_name() if inp.get_name() else "molecule"
        self._box = inp.get_box() if inp.get_box() else [0 for x in range(self._dim)]

        # Create output folder
        utils.mkdirp(link)


    ###############
    # Antechamber #
    ###############
    def job(self, name="", master=""):
        """Create job file to run with Antechamber. A shell file and a tleap
        file are created with all necessary commands to create a topology from a
        pdb file. For the conversion to GROMACS file format the python package
        **ParmED** needs to be installed.

        Parameters
        ----------
        name : string, optional
            Filename
        master : string, optional
            Set name for master job file if needed (practical for running
            multiple topology generations)
        """
        # Initialize
        link = self._link
        mol_name = self._name.lower()
        short = self._inp.get_short()
        name = name if name else self._name

        # Template directory
        package_dir = os.path.split(__file__)[0]+"/"

        # Job file
        file_in = package_dir+"templates/antechamber.job"
        file_out = link+name+".job"
        shutil.copy(file_in, file_out)

        utils.replace(file_out, "MOLNAME", mol_name)

        # Tleap file
        file_in = package_dir+"templates/antechamber.tleap"
        file_out = link+name+".tleap"
        shutil.copy(file_in, file_out)

        utils.replace(file_out, "MOLSHORTLOWER", mol_name)
        utils.replace(file_out, "MOLSHORT", short)
        utils.replace(file_out, "MOLNAME", mol_name)

        # Add to master run
        if master:
            fileMaster = open(link+master, "a")
            fileMaster.write("cd "+mol_name+" #\n")
            fileMaster.write("sh "+mol_name+".job #\n")
            fileMaster.write("cd .. #\n")
            fileMaster.write("echo \"Finished "+mol_name+"...\"\n")
            fileMaster.write("#\n")
            fileMaster.close()


    #############
    # Structure #
    #############
    def obj(self, name=""):
        """Save the molecule object or the dictionary of those using pickle.

        Parameters
        ----------
        name : string, optional
            Filename
        """
        # Initialize
        link = self._link
        link += name if name else self._name+".obj"

        # Save object
        utils.save(self._inp, link)

    def pdb(self, name="", use_atom_names=False):
        """Generate the structure file for the defined molecule in the **PDB**
        format.

        Parameters
        ----------
        name : string, optional
            Filename
        use_atom_names : bool, optional
            True to use atom names if they are defined, False to enumerate based
            on type
        """
        # Initialize
        link = self._link
        link += name if name else self._name+".pdb"

        # Open file
        with open(link, "w") as file_out:
            # Set counter
            num_a = 1
            num_m = 1

            # Run through molecules
            for mol in self._mols:
                atom_types = {}
                temp_res_id = 0
                # Run through atoms
                for atom in mol.get_atom_list():
                    # Process residue index
                    if not atom.get_residue() == temp_res_id:
                        num_m = num_m+1 if num_m < 9999 else 1
                        temp_res_id = atom.get_residue()

                    # Get atom type
                    atom_type = atom.get_atom_type()

                    # Create type dictionary
                    if atom_type not in atom_types:
                        atom_types[atom_type] = 1

                    # Set atom name
                    if use_atom_names and atom.get_name():
                        atom_name = atom.get_name()
                    else:
                        atom_name = atom_type+str(atom_types[atom_type])

                    # Write file
                    out_string = "HETATM"                  #  1- 6 (6)    Record name
                    out_string += "%5i" % num_a            #  7-11 (5)    Atom serial number
                    out_string += " "                      # 12    (1)    -
                    out_string += "%4s" % atom_name        # 13-16 (4)    Atom name
                    out_string += " "                      # 17    (1)    Alternate location indicator
                    out_string += "%3s" % mol.get_short()  # 18-20 (3)    Residue name
                    out_string += " "                      # 21    (1)    -
                    out_string += "%1s" % "A"              # 22    (1)    Chain identifier
                    out_string += "%4i" % num_m            # 23-26 (4)    Residue sequence number
                    out_string += " "                      # 27    (1)    Code for insertion of residues
                    out_string += "   "                    # 28-30 (3)    -
                    for i in range(self._dim):             # 31-54 (3*8)  Coordinates
                        out_string += "%8.3f" % (atom.get_pos()[i]*10)
                    out_string += "%6.2f" % 1              # 55-60 (6)    Occupancy
                    out_string += "%6.2f" % 0              # 61-66 (6)    Temperature factor
                    out_string += "           "            # 67-77 (11)   -
                    out_string += "%2s" % atom_type        # 78-79 (2)    Element symbol
                    out_string += "  "                     # 80-81 (2)    Charge on the atom

                    file_out.write(out_string+"\n")

                    # Process counter
                    num_a = num_a+1 if num_a < 99999 else 1
                    atom_types[atom_type] = atom_types[atom_type]+1 if atom_types[atom_type] < 99 else 1
                num_m = num_m+1 if num_m < 9999 else 1

            # End statement
            file_out.write("TER\nEND\n")

    def gro(self, name="", use_atom_names=False):
        """Generate the structure file for the defined molecule in the **GRO**
        format.

        Parameters
        ----------
        name : string, optional
            Filename
        use_atom_names : bool, optional
            True to use atom names if they are defined, False to enumerate based
            on type
        """
        # Initialize
        link = self._link
        link += name if name else self._name+".gro"

        # Open file
        with open(link, "w") as file_out:
            # Set title
            file_out.write("Molecule generated using the PoreMS package\n")

            # Number of atoms
            file_out.write("%i" % sum([x.get_num() for x in self._mols])+"\n")

            # Set counter
            num_a = 1
            num_m = 1

            # Run through molecules
            for mol in self._mols:
                atom_types = {}
                temp_res_id = 0
                # Run through atoms
                for atom in mol.get_atom_list():
                    # Process residue index
                    if not atom.get_residue() == temp_res_id:
                        num_m = num_m+1 if num_m < 99999 else 0
                        temp_res_id = atom.get_residue()

                    # Get atom type
                    atom_type = atom.get_atom_type()

                    # Create type dictionary
                    if atom_type not in atom_types:
                        atom_types[atom_type] = 1

                    # Set atom name
                    if use_atom_names and atom.get_name():
                        atom_name = atom.get_name()
                    else:
                        atom_name = atom_type+str(atom_types[atom_type])

                    # Write file
                    out_string = "%5i" % num_m              #  1- 5 (5)    Residue number
                    out_string += "%-5s" % mol.get_short()  #  6-10 (5)    Residue short name
                    out_string += "%5s" % atom_name         # 11-15 (5)    Atom name
                    out_string += "%5i" % num_a             # 16-20 (5)    Atom number
                    for i in range(self._dim):                    # 21-44 (3*8)  Coordinates
                        out_string += "%8.3f" % atom.get_pos()[i]

                    file_out.write(out_string+"\n")

                    # Process counter
                    num_a = num_a+1 if num_a < 99999 else 0
                    atom_types[atom_type] = atom_types[atom_type]+1 if atom_types[atom_type] < 999 else 0
                num_m = num_m+1 if num_m < 99999 else 0

            # Box
            out_string = ""
            for i in range(self._dim):
                out_string += "%.3f" % self._box[i]
                out_string += " " if i < self._dim-1 else "\n"

            file_out.write(out_string)

    def xyz(self, name="", use_atom_names=False):
        """Generate the structure file for the defined molecule in the **XYZ**
        format for running qm-simulations.

        Parameters
        ----------
        name : string, optional
            Filename
        use_atom_names : bool, optional
            True to use atom names if they are defined, False to enumerate based
            on type
        """
        # Initialize
        link = self._link
        link += name if name else self._name+".xyz"

        # Open output file and set title
        with open(link, "w") as file_out:
            # Header
            file_out.write("%i" % sum([x.get_num() for x in self._mols])+"\n"+"Energy = \n")

            # Run through molecules
            for mol in self._mols:
                # Run through atoms
                for atom in mol.get_atom_list():
                    # Write file
                    out_string = "%-2s" % atom.get_atom_type()  # 1- 2 (2)     Atom name
                    for i in range(self._dim):                  # 3-41 (3*13)  Coordinates
                        out_string += "%13.7f" % (atom.get_pos()[i]*10)

                    file_out.write(out_string+"\n")

    def lmp(self, name=""):
        """Generate the structure file for the defined molecule in the LAMMPS
        format. Assuming real units are used, the coordinates are in Angstroms.

        Parameters
        ----------
        name : string, optional
            Filename
        """
        # Initialize
        link = self._link
        link += name if name else self._name+".lmp"

        # Atom types
        atom_types = list(set(sum([[x.get_atom_type(i) for i in range(x.get_num())] for x in self._mols], [])))

        # Open file
        with open(link, "w") as file_out:
            # Set title
            file_out.write("Molecule generated using the PoreMS package\n\n")

            # Porperties section
            file_out.write("%i" % sum([x.get_num() for x in self._mols])+" atoms\n")
            file_out.write("%i" % len(atom_types)+" atom types\n")
            file_out.write("\n")

            # Box size - Periodic boundary conditions
            file_out.write("0.000 "+"%.3f" % (self._box[0]*10)+" xlo xhi\n")
            file_out.write("0.000 "+"%.3f" % (self._box[1]*10)+" ylo yhi\n")
            file_out.write("0.000 "+"%.3f" % (self._box[2]*10)+" zlo zhi\n")
            file_out.write("\n")

            # Masses
            file_out.write("Masses\n\n")
            for i, at in enumerate(atom_types):
                file_out.write("%i"%(i+1)+" "+"%8.3f"%db.get_mass(at)+"\n")
            file_out.write("\n")

            # Atoms
            file_out.write("Atoms\n\n")

            # Set counter
            num_a = 1
            num_m = 1

            # Run through molecules
            for mol in self._mols:
                temp_res_id = 0
                # Run through atoms
                for atom in mol.get_atom_list():
                    # Process residue index
                    if not atom.get_residue() == temp_res_id:
                        num_m = num_m+1
                        temp_res_id = atom.get_residue()

                    # Get atom type
                    atom_type_id = atom_types.index(atom.get_atom_type())+1

                    # Write atom line
                    out_string  = "%5i" % num_a + " "        #  Atom number
                    out_string += "%5i" % num_m + " "        #  Residue number
                    out_string += "%3i" % atom_type_id + " " #  Atom type
                    out_string += "%5i" % 0 + " "            #  Charge
                    for i in range(self._dim):               #  Coordinates
                        out_string += "%8.3f" % (atom.get_pos()[i]*10)
                        out_string += " " if i<self._dim-1 else ""
                    file_out.write(out_string+"\n")

                    # Process counter
                    num_a = num_a+1
                num_m = num_m+1


    ############
    # Topology #
    ############
    def top(self, name=""):
        """Store the **topology** file for a pore. A top file is created
        containing the itp-include for all molecules and the count of the
        different groups of the pore.

        Parameters
        ----------
        name : string, optional
            Filename
        """
        # Check input type
        if not isinstance(self._inp, Pore):
            print("Store: Unsupported input type for topology creation...")
            return

        # Initialize
        link = self._link
        link += name if name else self._name+".top"

        # Copy master topology file
        utils.copy(os.path.split(__file__)[0]+"/templates/topol.top", link)

        # Open file
        with open(link, "a") as file_out:
            # Include topology
            for mol_short in self._short_list:
                if mol_short not in ["SI", "OM", "SL", "SLG", "SLX"]:
                    file_out.write("#include \""+mol_short+".itp\"\n")

            file_out.write("\n")
            file_out.write("[ system ]\n")
            file_out.write("Pore-System Generated by the PoreMS Package\n")

            file_out.write("\n")
            file_out.write("[ molecules ]\n")

            # Number of atoms
            for mol_short in self._short_list:
                file_out.write(mol_short+" "+str(len(self._inp.get_mol_dict()[mol_short]))+"\n")

    def grid(self, name="", charges={"si": 1.28, "om": -0.64}):
        """Store the **grid.itp** file containing the necessary parameters and
        charges of the grid molecules.

        Parameters
        ----------
        name : string, optional
            Filename
        charges : dictionary, optional
            Dictionary of charges for silicon and oxygen atoms
        """
        # Initialize
        link = self._link
        link += name if name else "grid.itp"

        # Copy grid file
        utils.copy(os.path.split(__file__)[0]+"/templates/grid.itp", link)

        # Replace charges
        utils.replace(link, "CHARGEO", "%8.6f" % charges["om"])
        utils.replace(link, "CHARGESI", "%8.6f" % charges["si"])
