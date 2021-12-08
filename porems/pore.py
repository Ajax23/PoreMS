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
    """This class prepares and converts a given molecule block into a pore
    system.

    Methods are provided for preparing the surface according to procedures
    described in literature and for functionalizing the surface with any desired
    molecule objects.

    Parameters
    ----------
    block : Molecule
        Block molecule to be prepared as a pore
    matrix : Matrix
        Matrix object containing all bond information of given block molecule
    """
    def __init__(self, block, matrix):
        # Initialize
        self._dim = 3

        self._name = ""
        self._box = []

        self._block = block
        self._matrix = matrix

        self._num_in_ex = 0
        self._oxygen_ex = []

        self._mol_dict = {"block": {}, "in": {}, "ex": {}}


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
        while 3 in Counter(sum([self._matrix.get_matrix()[atom] for atom in self._matrix.bound(1)], [])).values():
            si_count = Counter(sum([self._matrix.get_matrix()[atom] for atom in self._matrix.bound(1)], []))
            for si, count in si_count.items():
                if count >= 3:
                    self._matrix.strip(si)

    def amorph(self, dist=0.05, accept=[0.1, 0.2], trials=100):
        """Make Structure amorphous.

        Parameters
        ----------
        dist : float, optional
            Maximal displacement distance
        accept : list, optional
            Acceptance criteria for allowed bond distances with partners
        trials : integer, optional
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

                # Displace if new bond length is in acceptance range
                if is_disp:
                    self._block.get_atom_list()[atom_id].set_pos(disp_pos)
                    break

    def exterior(self, gap):
        """Create an exterior surface for reservoir attachement. Using the
        periodic boundary distance as the indication, the bond between Si and O
        is broken and the silicon atom is saturated with and additional O atom
        using the original bond vector for orientation. At the end, the system
        is shifted by the given gap vector which represents the distance between
        bonds.

        Parameters
        ----------
        gap : list
            Gap between box boundary and crystal structure
        """
        # Initialize
        box = self._block.get_box()
        atom_list = self._block.get_atom_list()
        bound_list = self._matrix.get_matrix()

        # Get list of all si atoms
        si_list = [atom_id for atom_id, atom in enumerate(self._block.get_atom_list()) if atom.get_atom_type()=="Si"]

        # Run through silicon atoms
        add_list = []
        break_list = []
        for si in si_list:
            # Run through bound oxygen atoms
            for o in bound_list[si]:
                # Calculate bond vector
                bond_vector = [atom_list[si].get_pos()[dim]-atom_list[o].get_pos()[dim] for dim in range(3)]

                # Check if z dimension of bond - after rotation in pattern class - goes over boundary
                if abs(bond_vector[2]) > box[2]/2:
                    # Get bond direction
                    r = -0.155 if bond_vector[2] < 0 else 0.155

                    # Get bond angle
                    theta = geometry.angle_azi(bond_vector, is_deg=True)
                    phi = geometry.angle_polar(bond_vector, is_deg=True)
                    theta = theta-180 if r < 0 else theta

                    # Add new oxygen atom
                    self._block.add("O", si, r=r, theta=theta, phi=phi)

                    # Add list for new atoms
                    add_list.append([si, self._block.get_num()-1])

                    # Add to break bond list
                    break_list.append([si, o])

        # Move box to zero
        self._block.zero()
        self._block.translate(gap)

        # Break periodic bonds
        for bond in break_list:
            self._matrix.split(bond[0], bond[1])

        # Add bonds to new oxygens
        for bond in add_list:
            self._matrix.add(bond[0], bond[1])

        # Add atoms to exterior oxygen list
        self.prepare()
        self._oxygen_ex = self._matrix.bound(1)

    def sites(self):
        """Create binding site dictionary of the format

        .. math::

            \\boldsymbol{B}=\\begin{Bmatrix}
                si_1:&
                \\begin{Bmatrix}
                    \\text{o}: \\begin{pmatrix}o_{1,1},\\dots o_{1,n}\\end{pmatrix}&
                    \\text{type}: \\text{in/ex}&
                    \\text{state}: \\text{True/False}
                \\end{Bmatrix}\\\\
                \\vdots&\\vdots\\\\
                si_m:&
                \\begin{Bmatrix}
                    \\text{o}: \\begin{pmatrix}o_{m,1},\\dots o_{m,n}\\end{pmatrix}&
                    \\text{type}: \\text{in/ex}&
                    \\text{state}: \\text{True/False}
                \\end{Bmatrix}
            \\end{Bmatrix}

        with entries

        * **o** - Unsaturated oxygen atoms bound to silicon atom
        * **type** - Exterior "ex" or interior "in" binding site
        * **state** - available - True, used - False
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
            is_in = False
            is_ex = False
            for o in data["o"]:
                if o in self._oxygen_ex:
                    is_ex = True
                else:
                    is_in = True

            data["type"] = "ex" if is_ex else "in"

            if is_in and is_ex:
                self._num_in_ex += 1

            # State
            data["state"] = True


    #######################
    # Molecule Attachment #
    #######################
    def attach(self, mol, mount, axis, sites, amount, normal, scale=1, trials=1000, pos_list=[], site_type="in", is_proxi=True, is_random=True, is_rotate=False):
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
        amount : int
            Number of molecules to attach
        normal : function
            Function that returns the normal vector of the surface for a given
            position
        scale : float, optional
            Circumference scaling around the molecule position
        trials : integer, optional
            Number of trials picking a random site
        pos_list : list, optional
            List of positions (cartesian) to find nearest available binding site for
        site_type : string, optional
            Site type - interior **in**, exterior **ex**
        is_proxi : bool, optional
            True to fill binding sites in proximity of filled binding site
        is_random : bool, optional
            True to randomly pick a binding site from given list
        is_rotate : bool, optional
            True to randomly rotate molecule around own axis

        Returns
        -------
        mol_list : list
            List of molecule objects that are attached on the surface
        """
        # Check site type input
        if not site_type in ["in", "ex"]:
            print("Pore - Wrong attachement site-type...")
            return

        # Rotate molecule towards z-axis
        mol_axis = mol.bond(*axis)
        mol.rotate(geometry.cross_product(mol_axis, [0, 0, 1]), geometry.angle(mol_axis, [0, 0, 1]))
        mol.zero()

        # Search for overlapping placements - Calculate diameter and add carbon VdW-raidus (Wiki)
        if is_proxi:
            mol_diam = (max(mol.get_box()[:2])+0.17)*scale
            si_atoms = [self._block.get_atom_list()[atom] for atom in sites]
            si_dice = Dice(Molecule(inp=si_atoms), mol_diam, True)
            si_proxi = si_dice.find_parallel(None, ["Si", "Si"], [-mol_diam, mol_diam])
            si_matrix = {x[0]: x[1] for x in si_proxi}

        # Run through number of binding sites to add
        mol_list = []
        for i in range(amount):
            si = None
            # Find nearest free site if a position list is given
            if pos_list:
                pos = pos_list[i]
                min_dist = 100000000
                for site in sites:
                    if self._sites[site]["state"]:
                        length = geometry.length(geometry.vector(self._block.pos(site), pos))
                        if length < min_dist:
                            si = site
                            min_dist = length
            # Randomly pick an available site
            elif is_random:
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

                # Add molecule to molecule list and global dictionary
                mol_list.append(mol_temp)
                if not mol_temp.get_short() in self._mol_dict[site_type]:
                    self._mol_dict[site_type][mol_temp.get_short()] = []
                self._mol_dict[site_type][mol_temp.get_short()].append(mol_temp)

                # Remove bonds of occupied binding site
                self._matrix.strip([si]+self._sites[si]["o"])

                # Recursively fill sites in proximity with silanol and geminal silanol
                if is_proxi:
                    proxi_list = [sites[x] for x in si_matrix[sites.index(si)]]
                    if len(proxi_list) > 0:
                        mol_list.extend(self.attach(generic.silanol(), 0, [0, 1], proxi_list, len(proxi_list), normal, site_type=site_type, is_proxi=False, is_random=False))

        return mol_list

    def siloxane(self, sites, amount, normal, slx_dist=[0.507-1e-2, 0.507+1e-2], trials=1000, site_type="in"):
        """Attach siloxane bridges on the surface similar to Krishna et al.
        (2009). Here silicon atoms of silanol groups wich are at least 0.31 nm
        near each other can be converted to siloxan bridges, by removing one
        oxygen atom of the silanol groups and moving the other at the center of
        the two.

        Parameters
        ----------
        sites : list
            List of silicon ids of which binding sites should be picked
        amount : int
            Number of molecules to attach
        normal : function
            Function that returns the normal vector of the surface for a given
            position
        slx_dist : list
            Silicon atom distance bounds to search for parters in proximity
            [lower, upper]
        trials : integer, optional
            Number of trials picking a random site
        site_type : string, optional
            Site type - interior **in**, exterior **ex**
        """
        # Check site type input
        if not site_type in ["in", "ex"]:
            print("Pore - Wrong attachement site-type...")
            return

        # Create siloxane molecule
        mol = Molecule("siloxane", "SLX")
        mol.add("O", [0, 0, 0], name="OM1")
        mol.add("O", 0, r=0.09, name="OM1")
        mount = 0
        axis = [0, 1]

        # Rotate molecule towards z-axis
        mol_axis = mol.bond(*axis)
        mol.rotate(geometry.cross_product(mol_axis, [0, 0, 1]), geometry.angle(mol_axis, [0, 0, 1]))
        mol.zero()

        # Search for silicon atoms near each other
        si_atoms = [self._block.get_atom_list()[atom] for atom in sites]
        si_dice = Dice(Molecule(inp=si_atoms), slx_dist[1], False)
        si_proxi = si_dice.find_parallel(None, ["Si", "Si"], slx_dist)
        si_matrix = {x[0]: x[1] for x in si_proxi}

        # Run through number of siloxan bridges to add
        bond_matrix = self._matrix.get_matrix()
        mol_list = []
        for i in range(amount):
            # Randomly pick an available site pair
            si = []
            for j in range(trials):
                si_rand = random.choice(sites)
                # Check if binding site in local si-si matrix and if it contains binding partners
                if sites.index(si_rand) in si_matrix and si_matrix[sites.index(si_rand)]:
                    # Choose first binding partner in list
                    si_rand_proxi = sites[si_matrix[sites.index(si_rand)][0]]
                    # Check if binding partner is in local si-si matrix
                    if sites.index(si_rand_proxi) in si_matrix:
                        # Check if unbound states
                        if self._sites[si_rand]["state"] and self._sites[si_rand_proxi]["state"]:
                            # Check if binding site silicon atoms are already connected with an oxygen
                            is_connected = False
                            for atom_o in bond_matrix[si_rand]:
                                if atom_o in bond_matrix[si_rand_proxi]:
                                    is_connected = True
                            # Add to siloxane list if not connected
                            if not is_connected:
                                si = [si_rand, si_rand_proxi]
                                # End trials if appended
                                break

            # Place molecule on surface
            if si and self._sites[si[0]]["state"] and self._sites[si[1]]["state"]:
                # Create a copy of the molecule
                mol_temp = copy.deepcopy(mol)

                # Calculate center position
                pos_vec_halve = [x/2 for x in geometry.vector(self._block.pos(si[0]), self._block.pos(si[1]))]
                center_pos = [pos_vec_halve[x]+self._block.pos(si[0])[x] for x in range(self._dim)]

                # Rotate molecule towards surface normal vector
                surf_axis = normal(center_pos)
                mol_temp.rotate(geometry.cross_product([0, 0, 1], surf_axis), -geometry.angle([0, 0, 1], surf_axis))

                # Move molecule to mounting position and remove temporary atom
                mol_temp.move(mount, center_pos)
                mol_temp.delete(0)

                # Add molecule to molecule list and global dictionary
                mol_list.append(mol_temp)
                if not mol_temp.get_short() in self._mol_dict[site_type]:
                    self._mol_dict[site_type][mol_temp.get_short()] = []
                self._mol_dict[site_type][mol_temp.get_short()].append(mol_temp)

                # Remove oxygen atom and if not geminal delete site
                for si_id in si:
                    self._matrix.strip(self._sites[si_id]["o"][0])
                    if len(self._sites[si_id]["o"])==2:
                        self._sites[si_id]["o"].pop(0)
                    else:
                        del self._sites[si_id]
                    del si_matrix[sites.index(si_id)]

        return mol_list

    def fill_sites(self, sites, normal, site_type):
        """Fill list of given sites that are empty with silanol and geminal
        silanol molecules, respectively.

        Parameters
        ----------
        sites : list
            List of silicon ids
        normal : function
            Function that returns the normal vector of the surface for a given
            position
        site_type : string, optional
            Site type - interior **in**, exterior **ex**

        Returns
        -------
        mol_list : list
            List of molecule objects that are attached on the surface
        """
        mol_list = self.attach(generic.silanol(), 0, [0, 1], sites, len(sites), normal, site_type=site_type, is_proxi=False, is_random=False)

        return mol_list


    ###############
    # Final Edits #
    ###############
    def objectify(self, atoms):
        """Create molecule objects of specified list of atoms.

        Parameters
        ----------
        atoms : list
            List of atom ids

        Returns
        -------
        mol_list : list
            List of molecule objects
        """
        # Initialize
        mol_list = []

        # Run through all remaining atoms with a bond or more
        for atom_id in atoms:
            # Get atom object
            atom = self._block.get_atom_list()[atom_id]

            # Create molecule object
            if atom.get_atom_type() == "O":
                mol = Molecule("om", "OM")
                mol.add("O", atom.get_pos(), name="OM1")
            elif atom.get_atom_type() == "Si":
                mol = Molecule("si", "SI")
                mol.add("Si", atom.get_pos(), name="SI1")

            # Add to molecule list and global dictionary
            mol_list.append(mol)
            if not mol.get_short() in self._mol_dict["block"]:
                self._mol_dict["block"][mol.get_short()] = []
            self._mol_dict["block"][mol.get_short()].append(mol)

        # Output
        return mol_list

    def reservoir(self, size):
        """Create a reservoir with given length on each side of pore system.

        Parameters
        ----------
        size : float
            Reservoir size in nm
        """
        # Convert molecule dict into list
        mol_list = sum([x for x in self.get_mol_dict().values()], [])

        # Get zero translation
        min_z = 1000000
        max_z = 0
        for mol in mol_list:
            col_z = mol.column_pos()[2]
            min_z = min(col_z) if min(col_z) < min_z else min_z
            max_z = max(col_z) if max(col_z) > max_z else max_z

        # Translate all molecules
        for mol in mol_list:
            mol.translate([0, 0, -min_z+size])

        # Set new box size
        box = self._block.get_box()
        self.set_box([box[0], box[1], max_z-min_z+2*size])


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
        mol_dict = {}
        for site_type in self._mol_dict.keys():
            for key, item in self._mol_dict[site_type].items():
                if not key in mol_dict.keys():
                    mol_dict[key] = []
                mol_dict[key].extend(item)

        return mol_dict

    def get_site_dict(self):
        """Return the molecule dictionary with site type differentiation.

        Returns
        -------
        site_dict : dictionary
            Dictionary of all molecules with sorted sites
        """
        return self._mol_dict

    def get_num_in_ex(self):
        """Return the number of geminal si-binding sites that have one OH on the
        interior surface and on on the exterior surface.

        Returns
        -------
        num_in_ex : integer
            Dictionary of all molecules with sorted sites
        """
        return self._num_in_ex
