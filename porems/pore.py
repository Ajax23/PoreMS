################################################################################
# Pore Class                                                                   #
#                                                                              #
"""Function for editing a pore surface."""
################################################################################


import copy
import random

from collections import Counter

import porems.geometry as geometry
import porems.generic as generic

from porems.dice import Dice
from porems.molecule import Molecule


class Pore():
    """Pore class

    Parameters
    ----------

    Examples
    --------

    """
    def __init__(self, block, matrix):
        # Initialize
        self._dim = 3

        self._name = ""
        self._box = []

        self._block = block
        self._matrix = matrix

        self._mol_dict = {}


    ###########
    # Surface #
    ###########
    def prepare(self):
        """For the purpose of ensuring chemical propriety, the carved out
        surface needs to be processed based on a set of rules as proposed by
        Coasne et al. (2008).

        First, all unsaturated silicon atoms are to be removed. Next, silicon
        atoms with three unsaturated oxygen bonds are eliminated. Finally, now
        unbound oxygen atoms must be deleted. The resulting surface has fully
        saturated silicon atoms with a maximum of two unsaturated oxygen bonds.

        These unsaturated oxygens will be used as binding sites to connect
        molecules that are to be placed on the surface.
        """
        # Remove unsaturated silicon atoms
        for atom in self._matrix.bound(4, "lt"):
            if self._block.get_atom_type(atom)=="Si":
                self._matrix.strip(atom)

        # Remove silicon atoms with three unsaturated oxygen atoms
        si_count = Counter(sum([self._matrix.get_matrix()[atom] for atom in self._matrix.bound(1)], []))
        for si, count in si_count.items():
            if count == 3:
                self._matrix.strip(si)

    def amorph(self, dist=0.05, accept=[0.1, 0.2], trials=100):
        """Make Structure amorphous.

        Parameters
        ----------
        dist : float, optional
            Maximal displacement distance
        accept : list, optional
            Acceptance criteria for allowed bond distances with partners
        trials : integer
            Allowed number of trials per atom
        """
        # Get connectivity matrix
        connect = self._matrix.get_matrix()

        # Run through atoms that have bond partners
        for atom_id in self._matrix.bound(0, "gt"):
            # Create testing atom
            atom_temp = copy.deepcopy(self._block.get_atom_list()[atom_id])
            atom_temp_pos = atom_temp.get_pos()[:]

            # Run through trials
            for i in range(trials):
                # Create random displacement vector
                disp_vec = [random.uniform(-dist, dist) for x in range(self._dim)]
                disp_pos = [atom_temp_pos[x]+disp_vec[x] for x in range(self._dim)]

                # Displace test atom
                atom_temp.set_pos(disp_pos)

                # Calculate new bond lengths
                is_disp = True
                for bond in connect[atom_id]:
                    bond_length = geometry.length(self._block.bond(atom_id, bond))
                    if bond_length < accept[0] or bond_length > accept[1]:
                        is_disp = False
                        break

                # Displace if new bond length are in acceptance range
                if is_disp:
                    self._block.get_atom_list()[atom_id].set_pos(disp_pos)
                    break

    def sites(self, exterior=[]):
        """Create binding site dictionary of the format

        .. math::

            \\boldsymbol{B}=\\begin{Bmatrix}
                si_1:&
                \\begin{Bmatrix}
                    \\text{"o"}: \\begin{pmatrix}o_{1,1},\\dots o_{1,n}\\end{pmatrix}&
                    \\text{"type"}: \\text{"in"/"ex"}&
                    \\text{"state"}: \\text{True/False}
                \\end{Bmatrix}\\\\
                \\vdots
            \\end{bmatrix}

        with entries

        * **o** - Unsaturated oxygen atoms bound to silicon atom
        * **type** - Exterior "ex" or interior "in" binding site
        * **state** - available - True, used - False

        Parameters
        ----------
        exterior : list, optional
            List of oxygen atom ids that are on the exterior surface
        """
        # Get list of surface oxygen atoms
        oxygen_list = self._matrix.bound(1)
        connect = self._matrix.get_matrix()

        # Create binding site dictionary
        self._sites = {}
        for o in oxygen_list:
            if not connect[o][0] in self._sites:
                self._sites[connect[o][0]] = {"o": []}
            self._sites[connect[o][0]]["o"].append(o)

        # Fill other information
        for si, data in self._sites.items():
            # Site type
            if exterior:
                for o in data["o"]:
                    if o in exterior:
                        site_type = "ex"
                        break
                    else:
                        site_type = "in"
            else:
                site_type = "in"
            data["type"] = site_type

            # State
            data["state"] = True


    #######################
    # Molecule Attachment #
    #######################
    def attach(self, mol, mount, axis, sites, amount, normal, scale=1, trials=1000, is_proxi=True, is_random=True, is_rotate=False):
        """Attach molecules on the surface.

        Parameters
        ----------
        mol : Molecule
            Molecule object to attach
        mount : integer
            Atom id of the molecule that is placed on the surface silicon atom
        axis : list
            List of two atom ids of the molecule that define the molecule axis
        sites : list
            List of silicon ids of which binding sites should be picked
        normal : function
            Function that returns the normal vector of the surface for a given
            position
        scale : float, optional
            Circumference scaling around the molecule position
        trials : integer, optional
            Number of trials picking a random site
        is_proxi : bool, optional
            True to fill binding sites in proximity of filled bingin site
        is_random : bool, optional
            True to randomly pick a binding site from given list
        is_rotate : bool, optional
            True to randomly rotate molecule around own axis

        Returns
        -------
        mol_list : list
            List of molecule objects that are attached on the surface
        """
        # Calculate molecule diameter and add carbon Van der Waals radius (Wiki)
        mol_axis = mol.bond(*axis)
        mol.rotate(geometry.cross_product(mol_axis, [0, 0, 1]), geometry.angle(mol_axis, [0, 0, 1]))
        mol.zero()

        # Search for possible overlapping placements
        if is_proxi:
            mol_diam = (max(mol.get_box()[:2])+0.17)*scale
            si_atoms = [self._block.get_atom_list()[atom] for atom in sites]
            si_dice = Dice(Molecule(inp=si_atoms), mol_diam, True)
            si_proxi = si_dice.find_parallel(None, ["Si", "Si"], 0, mol_diam)
            si_matrix = {x[0]: x[1] for x in si_proxi}

        # Run through number of binding sites to add
        mol_list = []
        for i in range(amount):
            # Randomly pick an available site
            si = None
            if is_random:
                for j in range(trials):
                    si_rand = random.choice(sites)
                    if self._sites[si_rand]["state"]:
                        si = si_rand
                        break
            # Or use next binding site in given list
            else:
                si = sites[i] if i<len(sites) else None

            # Place molecule on surface
            if si is not None and self._sites[si]["state"]:
                # Disable binding site
                self._sites[si]["state"] = False

                # Create a copy of the molecule
                mol_temp = copy.deepcopy(mol)

                # Check if geminal
                if len(self._sites[si]["o"])==2:
                    mol_temp.add("O", mount, r=0.164, theta=45)
                    mol_temp.add("H", mol_temp.get_num()-1, r=0.098)
                    mol_temp.set_name(mol.get_name()+"g")
                    mol_temp.set_short(mol.get_short()+"G")

                # Rotate molecule towards surface normal vector
                surf_axis = normal(self._block.pos(si))
                mol_temp.rotate(geometry.cross_product([0, 0, 1], surf_axis), -geometry.angle([0, 0, 1], surf_axis))

                # Move molecule to mounting position
                mol_temp.move(mount, self._block.pos(si))

                # Add molecule to molecule list and dict
                mol_list.append(mol_temp)
                if not mol_temp.get_short() in self._mol_dict:
                    self._mol_dict[mol_temp.get_short()] = []
                self._mol_dict[mol_temp.get_short()].append(mol_temp)

                # Remove bonds of occupied binding site
                self._matrix.strip([si]+self._sites[si]["o"])

                # Recursively fill sites in proximity with silanol and geminal silanol
                if is_proxi:
                    proxi_list = [sites[x] for x in si_matrix[sites.index(si)]]
                    if len(proxi_list) > 0:
                        mol_list.extend(self.attach(generic.silanol(), 0, [0, 1], proxi_list, len(proxi_list), normal, is_proxi=False, is_random=False))

        return mol_list

    def attach_special(self):
        """Special attachment of molecules on the surface.

        TODO: Finish function.
        """
        return

    def attach_siloxane(self):
        """Attach siloxane bridges on the surface

        TODO: Finish function. dice search si with 0.27nm like in Krishna
        (2009), center oxygen atom, disable binding site, if geminal only remove
        used o entry from si site -> not geminal anymorem, matrix.strip()
        """
        return


    ###############
    # Final Edits #
    ###############
    def fill_sites(self, sites, normal):
        """Fill list of given sites that are empty with silanol and geminal
        silanol molecules respectively. Once

        Parameters
        ----------
        sites : list
            List of silicon ids
        normal : function
            Function that returns the normal vector of the surface for a given
            position

        Returns
        -------
        mol_list : list
            List of molecule objects that are attached on the surface
        """
        mol_list = self.attach(generic.silanol(), 0, [0, 1], sites, len(sites), normal, is_proxi=False, is_random=False)

        return mol_list

    def reservoir(self):
        """Create reservoir and center box.

        TODO: Finish function
        """
        return

    def objectify(self):
        """Create molecule objects of remaining grid atoms.

        Returns
        -------
        mol_list : list
            List of molecule objects
        """
        # Initialize
        mol_list = []

        # Run through all remaining atoms with a bond or more
        for atom_id in self._matrix.bound(0, "gt"):
            # Get atom object
            atom = self._block.get_atom_list()[atom_id]

            # Create molecule object
            if atom.get_atom_type() == "O":
                mol = Molecule("om", "OM")
                mol.add("O", atom.get_pos(), name="OM1")
            elif atom.get_atom_type() == "Si":
                mol = Molecule("si", "SI")
                mol.add("Si", atom.get_pos(), name="SI1")
            else:
                print("Pore: Unknown atom type...")
                return

            # Add to molecule list and global dictionary
            mol_list.append(mol)
            if not mol.get_short() in self._mol_dict:
                self._mol_dict[mol.get_short()] = []
            self._mol_dict[mol.get_short()].append(mol)

        # Output
        return mol_list


    ###########
    # Analyze #
    ###########
    def properties(self):
        """Calculate properties.

        TODO: As own class with table representation? - Takes volume and surface
        as inputs -> from shape class
        """
        return


    ##################
    # Setter Methods #
    ##################
    def set_name(self, name):
        """Set the pore name.

        Parameters
        ----------
        name : string
            Pore name
        """
        self._name = name

    def set_box(self, box):
        """Set the box size.

        Parameters
        ----------
        box : list
            Box size in all dimensions
        """
        self._box = box


    ##################
    # Getter Methods #
    ##################
    def get_name(self):
        """Return the pore name.

        Returns
        -------
        name : string
            Pore name
        """
        return self._name

    def get_box(self):
        """Return the box size of the pore.

        Returns
        -------
        box : list
            Box size in all dimensions
        """
        return self._box

    def get_block(self):
        """Return the block molecule.

        Returns
        -------
        block : Molecule
            Block molecule object
        """
        return self._block

    def get_sites(self):
        """Return the binding sites dictionary.

        Returns
        -------
        sites : dictionary
            Binding sites dictionary
        """
        return self._sites

    def get_mol_dict(self):
        """Return the dictionary of all molecules.

        Returns
        -------
        mol_dict : dictionary
            Dictionary of all molecules
        """
        return self._mol_dict
