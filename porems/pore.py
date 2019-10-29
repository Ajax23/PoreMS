################################################################################
# Pore Class                                                                   #
#                                                                              #
"""Extension of the molecule class for pores."""
################################################################################


import math
import copy
import random
import multiprocessing as mp

import porems.utils as utils

from porems.molecule import Molecule
from porems.verlet import Verlet
from porems.bonding import Bonding


class Pore(Molecule):
    """This class is an extension of the molecule class
    :class:`porems.molecule.Molecule` and the core class for creating
    silica pores.

    At the beginning a silica-oxygene grid is genreated from the smallest
    possible grid unit by duplication. Using the verlet list
    :class:`porems.verlet.Verlet` and bonding classes
    :class:`porems.bonding.Bonding`, this silica grid is prepared rapidly
    and binnding sites are provided.

    Here methods are provided for adding molecule groups to the
    surface of the generated pore. Furhermore a method is added to fill all
    unused binding sites with silanol and siloxan groups.

    All bound molecules are added to a global molecule list, in order to preserve
    individual molecule properties.

    Parameters
    ----------
    size : list, float
        Size of the silica-oxygene-grid, can be a float, if same size in all dimensions
    diam : float
        Pore diameter
    drill : str
        Axis for the pore drilling
    res : float
        Value to extend and translate the pore by on the drill side creating reservoirs
    vs : float
        Verlet-list box size
    is_pbc : bool
        True if periodic boundary conditions are needed
    is_time : bool
        True to print the used time for the single steps
    """
    def __init__(self, size, diam, drill, res=1, vs=0.4, is_pbc=True, is_time=False):
        # Call super class
        super(Pore, self).__init__()

        # Initialize
        self._dim = 3
        self._size = size if isinstance(size, list) else [size for s in range(self._dim)]
        self._diam = diam-{"x": 0.55, "y": 0.5, "z": 0.5}[drill]
        self._is_pbc = is_pbc
        self._is_time = is_time
        self._res = res
        self._vs = vs
        self._drill = drill

        # Define bond lengths
        self._si_grid = 0.1550
        self._si_oh = 0.1640
        self._oh = 0.0980

        # Set proximity range
        self._o_grid = 0.300

        # Define silica block information
        self._repeat = [0.506, 0.877, 1.240]
        self._gap = [0.126, 0.073, 0.155]
        self._num_rep = [round(self._size[i]/self._repeat[i]) for i in range(self._dim)]
        self._size = [self._repeat[i]*self._num_rep[i]+self._gap[i] for i in range(self._dim)]

        # Define charges
        self._q_si = 1.28
        self._q_o = -0.64
        self._q_oh = -0.74+0.42

        # Define lists
        self._site = []
        self._grid = []

        # Time management
        self._t_tot = 0

        # Create silica pore
        t = utils.tic()
        self._block(0, self._foundation())              # Build block
        self._orientation()                             # Rotate drill axis
        self.translate(self._gap)                       # Translate gap
        self._center = self._focal()                    # Find focal point
        self._t_tot += utils.toc(t, "Build  ", is_time)

        # Verlet and bonding
        t = utils.tic()
        self._verlet = Verlet(self, vs, is_pbc)         # Create verlet boxes
        self._bonding = Bonding(self._verlet)           # Create bond matrix
        self._t_tot += utils.toc(t, "Matrix ", is_time)

        # Prepare pore
        t = utils.tic()
        self._box = self.get_box()                     # Set box
        self._bonding.attach()                         # Preperare sides
        self._bonding.drill(self._center, self._diam)  # Drill pore
        self._bonding.prepare()                        # Prepare pore surface
        self._t_tot += utils.toc(t, "Prepare", is_time)

        # Find binding sites
        t = utils.tic()
        self.zero()
        self._box = self.get_box()                     # Reset box size
        self._bind()                                   # Create site array
        self._proxi()                                  # Find sites in proximity
        self._diam = self._diameter()                  # Calculate new diameter
        self._t_tot += utils.toc(t, "Binding", is_time)


    ##########################
    # Private Methods - Pore #
    ##########################
    def _foundation(self):
        """In this method the smallest possible grid unit is created in a way,
        that duplication is simply done by copying the molecule without the need
        of deleting atoms in the end.

        Returns
        -------
        block : Molecule
            Smallest possible pore unit
        """
        # Initialize
        dim = self._dim
        dist = self._si_grid
        angle = 2*math.atan(math.sqrt(2))*180/math.pi

        # Create pore
        h00 = Molecule()
        h00.add("Si", [0, 0, 0])
        h00.add("O", 0, r=dist, theta=angle, phi=60)
        h00.add("Si", 1, bond=[0, 1], r=dist)
        h00.add("O", 2, r=dist, theta=180,  phi=0)
        h00.add("Si", 3, bond=[2, 3], r=dist)
        h00.add("O", 4, r=dist, theta=angle, phi=300)
        h00.add("Si", 5, bond=[4, 5], r=dist)
        h00.add("O", 6, r=-dist, theta=angle, phi=60)
        h00.add("Si", 7, bond=[6, 7], r=dist)
        h00.add("O", 8, r=-dist, theta=180,  phi=0)
        h00.add("Si", 9, bond=[8, 9], r=dist)
        h00.add("O", 10, r=-dist, theta=angle, phi=300)

        h00.rotate("y", 90)
        h00.rotate("z", 90)
        h00.rotate("y", 180)
        h00.rotate("x", h00.angle(h00.bond(8, 6)[0], self._axis("y")))
        h00.rotate("x", -90)

        h01 = copy.deepcopy(h00)
        h01.add("O", 2, r=dist)
        h01.add("O", 6, r=dist)
        h01.add("O", 10, r=dist)

        h02 = copy.deepcopy(h00)
        h02.move(0, h01.pos(12))
        h02.translate([0, 0, dist])
        h02.delete([1, 2, 3, 4, 5, 6, 7])

        h03 = copy.deepcopy(h00)
        h03.move(0, h01.pos(14))
        h03.translate([0, 0, dist])
        h03.delete([2, 3, 4, 5, 6, 7, 8, 9, 10, 11])

        hexa = Molecule(inp=[h01, h02, h03])

        # Calculate xrepeat
        z = [hexa]
        for i in range(10):
            z.append(copy.deepcopy(hexa))

        z[1].rotate("x", 180)
        z[1].move(18, z[0].pos(18))
        z[1].translate([0, 0, dist*2])
        z[1].add("O", z[0].pos(18), r=dist)

        z[2].rotate("x", 180)
        z[2].move(0, z[0].pos(18))

        z[3].move(4, z[0].pos(15))

        z[4].rotate("x", 180)
        z[4].move(16, z[0].pos(18))
        z[4].delete([21, 19, 15, 20, 14, 12, 10, 2, 11, 1, 0])

        z[5].move(6, z[1].pos(16))
        z[5].delete([21, 19, 15, 20, 14, 12, 10, 2, 11, 1, 0])

        # Create repeat block
        block = Molecule(inp=z)
        block.overlap()

        # Determine number of needed blocks
        block.zero()

        # Delete overlapping atoms
        block.delete([4, 3, 2, 12, 15, 45, 46, 55, 57, 63, 26, 25, 24, 34, 37, 61])  # x-Axis
        block.delete([50, 48, 36, 40, 41])                               # y-Axis
        block.delete([45, 20, 19, 21, 22, 23, 24, 25, 17, 18])                # z-Axis

        # Move to zero
        block.zero()

        return block

    def _block(self, dim, block):
        """Reculsively duplicate and translate a given molecule block in all
        given dimensions. The duplication stops if the next added block would
        create the pore longer than specified in the constructor.

        Parameters
        ----------
        dim : int
            Repeat dimensions
        block : Molecule
            Molecule unit to be duplicated
        """
        if dim < self._dim:
            p = []

            for i in range(self._num_rep[dim]):
                temp = copy.deepcopy(block)
                vec = [(i+1)*self._repeat[dim] if j == dim else 0 for j in range(self._dim)]

                temp.translate(vec)

                p.append(temp)

            self._block(dim+1, Molecule(inp=p))
        else:
            self._append(block)

    def _orientation(self):
        """Rotate pore orientation, so that the specified drill axis becomes
        the z-axis.
        """
        # Initialize
        drill = self._drill
        gap = self._gap
        size = self._size

        # Rotate pore
        if drill == "x":
            self.rotate("y", 90)
        elif drill == "y":
            self.rotate("x", 90)

        # Set zero
        self.zero()

        # Update gap and size lists
        if drill == "x":
            self._gap = [gap[2], gap[1], gap[0]]
            self._size = [size[2], size[1], size[0]]
        elif drill == "y":
            self._gap = [gap[0], gap[2], gap[1]]
            self._size = [size[0], size[2], size[1]]


    ###################################
    # Private Methods - Binding sites #
    ###################################
    def _bind(self):
        """Create binding site matrix. This matrix :math:`\\boldsymbol{B}`
        is mostly generated by the bonding class :class:`poresms.bonding.Bonding`

        .. math::

            \\boldsymbol{B}=\\begin{bmatrix}
                o_0&s_0&p_0&u_0&t_0&g_0&w_0\\\\
                o_1&s_1&p_1&u_1&t_1&g_1&w_1\\\\
                \\vdots\\\\
                o_b&s_b&p_b&u_b&t_0&g_b&w_b\\\\
            \\end{bmatrix}

        with binding site number :math:`b` and entries

        0. Oxygen atom id :math:`o_{k=0,\\dots,b}`
        1. Silica atom id :math:`s_k`
        2. List of pointers :math:`p_k` to binding sites that are in proximity
        3. State :math:`u_k`: avilable - 0, used - 1, is in proximity - 2
        4. Type :math:`t_k`: inside the pore - 0, outsite the pore - 1
        5. Pointer :math:`g_k` to geminal binding site
        6. Molecule position :math:`w_k` in the global write list
        """
        # Initialize
        site = []

        # Get data
        osi = self._bonding.site()

        # Create bonding site list
        for os in osi:
            site.append([os[0], os[1], [], 0, os[2], os[3], -1])

        self._site.extend(site)

    def _proxi(self):
        """Using verlet lists, search for binding sites that are in proximity,
        by looking for oxygene atoms that are near each other.

        Fill the result in the binding site matrix in the proximity entry.
        """
        # Initialize
        site = self._site

        # Create local verlet list
        atoms = [o[0] for o in site]
        verlet = Verlet(self._temp(atoms), self._vs, self._is_pbc)

        # Find oxygenes in proximity
        box_list = [i for i in range(len(verlet.get_box()[1]))]
        oo = verlet.find_parallel(box_list, ["O", "O"], self._o_grid, 10e-2)
        ooCol = utils.column(oo)

        # Append proxi list to sites
        for i in range(len(site)):
            os = site[i]
            o = oo[ooCol[0].index(i)][1] if i in ooCol[0] else []

            os[2].extend(o)

    def _random(self, site_type, rate, inp, counter, is_proxi=False):
        """Get a random binding site list for a given percentage of the available binding sites.

        If a site is available (state is 0), it is added to the output-list and
        the state entry for this site is set to 1.
        Then the states of all binding sites in the proximity list are set to 2.

        In case the random site is not available, this run is repeated. A break-counter
        prevents the occurance of an endless loop.

        If *is_proxi* is set to True, then binding sites nearest to each other
        are determined. the return value then consists of lists containing the
        ids of both sites. Geminal pairs are not allowed, since creating the
        topology is too complex. This however might be added in the future.

        Parameters
        ----------
        site_type : int
            Type of the sites inside-0, outside-1
        rate : float
            Rate of how many binding sites
        inp : str
            Input type of the rate
        counter : int
            Number of attempts before breaking the randoming loop
        is_proxi : bool
            True to search for site pairs in proximity

        Returns
        -------
        rand : list
            List of pointers to random binding sites.
        """
        # Initialize
        site = self._site
        temp = utils.column(site)
        site_len = temp[4].count(site_type)
        id_list = [i for i in range(len(site))]

        # Calculate random element number
        if inp == "num":
            ran_num = rate
        elif inp == "percent":
            ran_num = site_len*rate
            ran_num /= 100
            ran_num = math.floor(ran_num)

        # Get list of all partners
        if is_proxi:
            site_min = []
            for i in range(len(site)):
                length = [self._length(self._vector(site[i][0], site[j][0])) if not i ==
                          j and site[i][4] == site[j][4] else 100000 for j in range(len(site))]

                site_min.append([length.index(min(length)), min(length)])

        # Get random elements from site list
        ran_list = []
        count = 0
        while len(ran_list) < ran_num and count < counter:
            # Get random site index
            rand = random.choice(id_list)

            # Check binding site
            is_app = (is_proxi and site[rand][4] == site_type and
                      site[rand][3] == 0 and site[site_min[rand][0]][3] == 0 and
                      site[rand][5] == None and site[site_min[rand][0]][5] == None)

            is_app = site[rand][4] == site_type and site[rand][3] == 0 if not is_proxi else is_app

            # Append to return list
            if is_app:
                # Reset counter
                count = 0

                # Append
                ran_list.append([rand, site_min[rand][0]] if is_proxi else rand)

                # Set site and proximity to used
                self._close([rand, site_min[rand][0]] if is_proxi else rand)
            else:
                # Add counter
                count += 1

        return ran_list

    def _close(self, sites):
        """This function cloese the specified site.

        Parameters
        ----------
        sites : int, list
            identifiers of sites to be closed
        """
        # Initialize
        site = self._site

        # Process input
        sites = [sites] if isinstance(sites, int) else sites

        # Close sites
        for x in sites:
            site[x][3] = 1
            for prox in site[x][2]:
                if not site[prox][3] == 1:
                    site[prox][3] = 2


    ##################################
    # Private Methods - Add Molecule #
    ##################################
    def _couple(self, mol, sio, site_type, rate, inp="percent", sites=None, counter=1000):
        """Add a molecule to the pore binding site.

        The *sio* input defines the silica atom in the first entry and the
        directional vector with the other two entries.

        First the molecule is rotated around its directional vector randomly to
        simulate a state closer to reality. Second the molecule is rotated so that
        its directional vector and the binding site vector match. One final rotation
        is applied for molecules on the inside, so that the directional vector
        points towards the center of the pore.

        The molecule is then moved, so that its silica atom is placed on the silica
        of the binding site. The oxygen is then placed ontop of the binding
        site oxygen, in order to keep the sites geometry.

        Finally the atoms of the binding sites are added to the list of atoms
        to be deleted at the end of the construction.

        These molecules are then sorted by geminal binding sites and then added
        to the global molecule list.

        :TODO: Make molecule move independent of sio vector - line 514

        Parameters
        ----------
        mol : Molecule
            Molecule to be added
        sio : list
            List containing the silica id and directional vector of the molecule
        site_type : int
            Type of the sites inside-0, outside-1
        rate : float
            Rate of how many binding sites
        inp : str
            Input type of the rate
        sites : int, list
            Listindex of binding sites
        counter : int
            Number of attempts before breaking the randoming loop
        """
        # Initialize
        dim = self._dim
        site = self._site

        # Process user input
        sites = [sites] if isinstance(sites, int) else sites

        # Get random list or given list
        ran_list = self._random(site_type, rate, inp, counter) if sites is None else sites

        # Molecule and center vector
        vec_m = mol.bond(sio[1], sio[2])[0]
        center = [self._center[i] for i in range(dim-1)]

        # Add molcule
        write_list = [[], []]
        for i in sorted(ran_list, reverse=True):
            # Initialize
            o = site[i]

            # Copy molecule
            temp = copy.deepcopy(mol)

            # Calculate osi vector
            vec_o = self._vector(o[1], o[0])
            vec_c = center[:]
            vec_c.append(self.pos(o[0])[dim-1])
            vec_c = self._vector(self.pos(o[1]), vec_c)

            # Random axis rotation
            temp.rotate(vec_m, random.randint(1, 180))

            # Adjust molecule vector to osi
            temp.rotate(self._cross(vec_m, vec_o), self.angle(vec_m, vec_o))

            # Adjust molecule towards center
            if site_type == 0:
                temp.rotate(self._cross(vec_c, vec_o), -self.angle(vec_c, vec_o))

            # Move molecule to position
            temp.move(sio[1], self.pos(o[0]))

            # Move silica
            temp.put(sio[0], self.pos(o[1]))

            # Check if geminal
            is_gem = not o[5] == None

            name = temp.get_name()+"g" if is_gem else temp.get_name()
            short = temp.get_short()+"g" if is_gem else temp.get_short()
            temp.set_name(name)
            temp.set_short(short)

            # Remove atoms
            self._bonding.remove([o[0], o[1]])

            # Add molecule to pore
            write_id = 1 if is_gem else 0
            write_list[write_id].append([temp, i])

        # Add to write
        for write in write_list:
            for i in range(len(write)):
                site[write[i][1]][6] = len(self._write)
                self._write.append(write[i][0])

    def _couple_proxi(self, mol, sio, orient, rate, inp="percent", counter=1000):
        """Add a molecule to two pore binding sites.

        The *sio* input defines the silica and oxygen atom pair in two seperate
        lists.

        First the molecule is rotated, so that the oxygens are ontop of the
        binding site oxygenes. Second the molecule is rotated so that the
        directional vector points to the center of the pore.

        The molecule is then moved, so that its oxygen atom is placed on the oxygen
        of the first binding site then translated, so the molecule is in the center
        of both sites.

        Finally the atoms of the binding sites are added to the list of atoms
        to be deleted at the end of the construction.

        Parameters
        ----------
        mol : Molecule
            Molecule to be added
        sio : list
            List containing two silica-oxygen vectors
        orient : list
            Molecule direction vector
        rate : float
            Rate of how many binding sites
        inp : str
            Input type of the rate
        counter : int
            Number of attempts before breaking the randoming loop
        """
        # Initialize
        dim = self._dim
        site = self._site

        # Get random list
        ran_list = self._random(0, rate, inp, counter, is_proxi=True)

        # Molecule and center vector
        center = [self._center[i] for i in range(dim-1)]

        # Add molcule
        for i in sorted(ran_list, reverse=True):
            # Initialize
            s = [site[i[0]], site[i[1]]]

            # Copy molecule
            temp = copy.deepcopy(mol)

            # Rotate towards binding site
            vec_mo = temp.bond(sio[0][1], sio[1][1])[0]  # Molecule oxygenes
            vec_so = self._vector(s[0][0], s[1][0])     # Binding site oxygenes
            vec_m = temp.bond(orient[0], orient[1])[0]  # Molecule orientation

            temp.rotate(vec_m, self.angle(vec_mo, vec_so))

            # Rotate towards central axis
            cent_o = [(self.pos(s[0][0])[x]+self.pos(s[1][0])[x])/2 for x in range(dim)]
            cent_p = center[:]+[cent_o[-1]]

            vec_m = temp.bond(orient[0], orient[1])[0]  # Molecule orientation
            vec_c = self._vector(cent_o, cent_p)       # Binding site center

            temp.rotate(self._cross(vec_m, vec_c), self.angle(vec_m, vec_c))

            # Put oxygenes ontop of each other
            temp.move(sio[0][1], self.pos(s[0][0]))

            # Move to center
            cent_m = [(temp.pos(sio[0][1])[x]+temp.pos(sio[1][1])[x])/2 for x in range(dim)]
            temp.translate(self._vector(cent_m, cent_o))

            # Move silica atoms
            temp.put(sio[0][0], self.pos(s[0][1]))
            temp.put(sio[1][0], self.pos(s[1][1]))

            # Remove atoms
            self._bonding.remove([s[0][0], s[0][1], s[1][0], s[1][1]])

            # Add to write
            self._write.append(temp)

    def _special(self, mol, sio, num, symmetry):
        """Add a molecule in a specified orientation with a specific amount.
        Currently following symmetry orientations to each other are available

        * **random** - No symmetry, random gaussian placement
        * **point** - Point symmetry
        * **mirror** - Mirror symmetry

        Parameters
        ----------
        mol : Molecule
            Molecule to be added
        sio : list
            List containing the silica id and directional vector of the molecule
        num : int
            Number of molecules to be added
        symmetry : str
            Molecule symmetry orientation to each other
        """
        # Initialize
        diam = self._diam
        site = self._site
        center = self._center

        # Process input
        if symmetry not in ["random", "point", "mirror"]:
            print("Symmetry type not supported...")
            return

        elif not symmetry == "random" and not num == 0:
            # Calculate placment points
            pos_list = []
            length = center[2]*2
            dist = length/num
            start = dist/2

            for i in range(num):
                if symmetry == "point":
                    coeff = -1 if i % 2 == 0 else 1
                elif symmetry == "mirror":
                    coeff = 1

                x = center[0]+coeff*diam/2
                y = center[1]
                z = start+dist*i

                pos_list.append([x, y, z])

            # Find nearest site
            site_min = [[0, None] for i in range(num)]
            for i in range(len(site)):
                length = [self._length(self._vector(self.pos(site[i][0]), pos)) for pos in pos_list]

                for j in range(num):
                    if site_min[j][1] is None or length[j] < site_min[j][1]:
                        site_min[j] = [i, length[j]]

            sites = [x[0] for x in site_min]

            self._close(sites)

        else:
            sites = None

        self._couple(mol, sio, 0, num, inp="num", sites=sites)


    ################################
    # Private Methods - Final Mols #
    ################################
    def _silanol(self, site_range=None):
        """Add hydrogene atoms to all open binding sites. Add rotate the direction
        of the binding site to point to the pore center,
        in case the site is on the inside. This is done in the same manner as
        function :func:`couple`.

        Parameters
        ----------
        site : list, None
            Site list to go trhough

        Returns
        -------
        write_list : list
            List of Molecule objects and ids.
        """
        # Initialize
        site = self._site
        dim = self._dim
        box = self._box
        center = [self._center[i] for i in range(dim-1)]
        write_list = [[], []]

        # Set siterange
        site_r = [0, len(site)] if site_range is None else site_range

        # Define temporary molecule object
        temp_mol = Molecule()

        # Change bond length and add hydrogen
        for i in range(site_r[0], site_r[1]):
            if not site[i][3] == 1:
                # Initialize
                o = site[i]
                o[3] = 1
                temp = copy.deepcopy(temp_mol)

                # Change SiO bond length
                if self.bond(o[0], o[1])[1] < box[2]/2:
                    self.part_move([o[1], o[0]], o[0], self._si_oh)

                # Inside pore
                if o[4] == 0:
                    temp.add("O", [0, 0, 0])
                    temp.add("H", 0, r=self._oh)

                    posC = center[:]
                    posC.append(self.pos(o[0])[dim-1])

                    vec_m = temp.bond(0, 1)[0]
                    vec_c = self._vector(self.pos(o[1]), posC)
                    vec_o = self._vector(o[1], o[0])

                    temp.rotate(self._cross(vec_m, vec_o), self.angle(vec_m, vec_o))
                    temp.rotate(self._cross(vec_c, vec_o), -self.angle(vec_c, vec_o))

                    temp.move(0, self.pos(o[0]))
                    temp.add("Si", self.pos(o[1]))

                    # Sort
                    temp_pos_o = temp.pos(0)
                    temp_pos_h = temp.pos(1)
                    temp.delete([0, 1])
                    temp.add("O", temp_pos_o)
                    temp.add("H", temp_pos_h)

                # Outside pore
                else:
                    angle = -180 if self.pos(o[0])[2] < box[2]/2 else 0

                    temp.add("Si", self.pos(o[1]))
                    temp.add("O", self.pos(o[0]))
                    temp.add("H", 1, r=self._oh, theta=angle)

                # Check if geminal
                is_gem = not o[5] == None

                # Add Charge
                temp.set_charge(self._q_si+self._q_oh)

                # Add to write
                sl_name = "slg" if is_gem else "sl"
                sl_short = "SLG" if is_gem else "SL"
                temp.set_name(sl_name)
                temp.set_short(sl_short)

                write_id = 1 if is_gem else 0
                write_list[write_id].append([temp, i])

        # Add to write and add site - molecule connection
        if site_range is None:
            for write in write_list:
                for w in write:
                    self._bonding.remove([site[w[1]][1], site[w[1]][0]])
                    site[w[1]][6] = len(self._write)
                    self._write.append(w[0])
        else:
            return write_list

    def _silanol_parallel(self):
        """Parallelized function :func:`_silanol`.
        """
        # Initialize
        np = mp.cpu_count()
        sites = self._site

        # Split site list
        site_len = math.floor(len(sites)/np)
        site_list = []
        for i in range(np):
            if i == np-1:
                site_list.append([site_len*i, len(sites)])
            else:
                site_list.append([site_len*i, site_len*(i+1)])

        # Paralellize
        pool = mp.Pool(processes=np)
        results = pool.map_async(self._silanol, site_list)
        pool.close()

        write_list = [[], []]
        for result in results.get():
            write_list[0].extend(result[0])
            write_list[1].extend(result[1])

        # Add to write and add site - molecule connection
        for write in write_list:
            for w in write:
                self._bonding.remove([sites[w[1]][1], sites[w[1]][0]])
                sites[w[1]][6] = len(self._write)
                self._write.append(w[0])

    # Transform two SiO to siloxan bridges
    def _siloxan(self, rate, inp="percent", counter=1000):
        """Add siloxan bridges to the pore by selecting two binding sites in
        proximity removing both oxygens and placing one at center of
        the two removed atoms.

        Parameters
        ----------
        rate : float
            Rate of binding sites to be editied
        inp : str
            Rate input type
        counter : int
            Number of attempts before breaking the randoming loop
        """
        # Initialize
        site = self._site
        dim = self._dim

        # Exit if rate is zero
        if rate == 0:
            return

        # Get random list
        ran_list = self._random(0, rate, inp, counter, is_proxi=True)

        # Molecule and center vector
        center = [self._center[i] for i in range(dim-1)]

        # Define temporary molecule object
        temp_mol = Molecule()

        # Add molcule
        for i in sorted(ran_list, reverse=True):
            # Initialize
            s = [site[i[0]], site[i[1]]]

            # Get center position
            cent_o = [(self.pos(s[0][0])[x]+self.pos(s[1][0])[x])/2 for x in range(dim)]
            cent_s = [(self.pos(s[0][1])[x]+self.pos(s[1][1])[x])/2 for x in range(dim)]

            # Create molecule
            temp = copy.deepcopy(temp_mol)
            temp.set_name("slx")
            temp.set_short("SLX")
            temp.add("OM", cent_o)
            temp.translate(self._vector(cent_s, cent_o))
            temp.add("O", self.pos(s[0][0]))
            temp.add("O", self.pos(s[1][0]))

            # Rotate towards central axis
            cent_p = center[:]+[cent_o[-1]]

            vec_m = self._vector(cent_o, temp.pos(0))
            vec_c = self._vector(cent_o, cent_p)

            temp.rotate(self._cross(vec_m, vec_c), self.angle(vec_m, vec_c))

            # Put oxygenes ontop of each other
            temp.move(1, self.pos(s[0][0]))

            # Remove temporary oxygen atoms
            temp.delete([1, 2])

            # Remove atoms
            self._bonding.remove([s[0][0], s[1][0]])

            # Add to write
            self._write.append(temp)


    #################################
    # Private Methods - Final Edits #
    #################################
    def _connect(self):
        """Concatenate two geminal molecules to one molecule.
        """
        # Initialize
        site = self._site
        write = self._write
        dim = self._dim
        gemi = []

        # Check for geminal sites
        for i in range(len(site)):
            if site[i][5] is not None:
                gem = [site[i][6], site[site[i][5]][6]]
                gem_r = [site[site[i][5]][6], site[i][6]]

                if not gem_r in gemi:
                    gemi.append(gem)

        # Append atoms
        pop = []
        for g in gemi:
            id_a = 0 if write[g[1]].get_name() == "slg" else 1
            id_b = abs(id_a-1)

            mol_a = write[g[id_a]]
            mol_b = write[g[id_b]]

            num = mol_b.get_num()

            for i in range(num):
                if not mol_b.get_type(i) == "Si":
                    mol_a.add(mol_b.get_type(i), mol_b.pos(i))

            pop.append(g[id_b])

            # Add charge
            write[g[id_a]].set_charge(write[g[id_a]].get_charge()+self._q_oh)

        # Delete from write list
        for d in sorted(pop, reverse=True):
            write.pop(d)

    def _objectify(self):
        """Move all remaining grid silica and oxygen atoms to individual molecules
        for writing the structure file and add these molecules to the global
        molecule list after sorting.
        """
        # Initialize
        data = self._data
        dim = self._dim
        write_list = [[], []]

        # Define temporary molecule object
        temp_mol = Molecule()

        # Create OM and SI mols
        for i in range(len(data[0])):
            temp = copy.deepcopy(temp_mol)

            # Check if oxygene
            is_oxy = self.get_type(i) == "O"

            # Set temporary vars
            atom = "OM" if is_oxy else "SI"
            name = "om" if is_oxy else "si"
            short = "OM" if is_oxy else "SI"
            charge = self._q_o if is_oxy else self._q_si
            write_id = 0 if is_oxy else 1

            # Testing
            if self.get_type(i) == "R":
                atom = "R"
                name = "test"

            # Creat unique object and add to write
            temp.add(atom, self.pos(i))
            temp.set_name(name)
            temp.set_short(short)
            temp.set_charge(charge)
            write_list[write_id].append(temp)

        # Add mols to write
        self._write.pop(0)

        for write in write_list:
            self._write = write + self._write

    def _sort(self):
        """Sort molecules in order to prevent multiple molecule appearances in
        the structure file.
        """
        # Initialize
        write = self._write
        write_new = []

        # Get unique mols
        unique_mols = []
        for mol in write:
            if mol.get_name() not in unique_mols:
                unique_mols.append(mol.get_name())

        # Sort molecules
        for molName in unique_mols:
            for mol in write:
                if mol.get_name() == molName:
                    write_new.append(mol)

        self._write = write_new

    def _position(self):
        """Tranlate the pore into position, and consider the periodic repeat gap.
        Hereby the coordinates are moved by halve the gap distance and box box is
        extended by the other halve. Additionally, to prevent molecule overlap and
        perodic movement, the drill side is extended and translated by a fixed value,
        thus creating reservoirs.
        """
        # Initialize
        dim = self._dim
        box = self._box
        gap = self._gap
        write = self._write

        # Get vector to coordinate zero
        vec = Molecule(inp=[x for x in write if x.get_name() in ["si", "om"]]).zero()

        # Position pore
        for x in write:
            x.translate([vec[i]+gap[i]/2 for i in range(dim)])
            x.translate([self._res if i == 2 else 0 for i in range(dim)])

        # Extend box
        coord = Molecule(inp=[x for x in write if x.get_name() in ["si", "om", "ox"]]).get_box()

        for i in range(dim):
            box[i] = coord[i] + gap[i]/2
        box[2] += self._res

    def _overlap(self):
        """Method for checking of silica or oxygene atoms are overlapping
        using verlet lists. Duplicate atoms will be printed.
        """
        # Initialize
        temp = Molecule()
        for write in self._write:
            temp._append(write)

        # Create verlet
        verlet = Verlet(temp, self._vs, True)

        # Check silica and oxygene
        si = verlet.find_parallel(None, ["Si", "Si"], 0, 10e-3)
        oo = verlet.find_parallel(None, ["O", "O"], 0, 10e-3)

        # Output
        for s in si:
            if len(s[1]) > 0:
                print(s)

        for o in oo:
            if len(o[1]) > 0:
                print(o)

    def _excess(self, is_rand=True):
        """Adjust grid atoms charges to compensate for excess charges on the pore.

        Parameters
        ----------
        is_rand : bool
            True to randomly distribute rounding error charge on block oxygen atoms
        """
        # Initialize
        mols = self._write
        grid = [x.get_name() for x in self._grid]+["sl"]
        grid += [x+"g" for x in grid]

        # Get excess charge
        charge = 0.0
        num_o = num_si = num_g = 0
        round_o = 6
        for mol in mols:
            charge += mol.get_charge()

            if mol.get_name() == "om":
                num_o += 1
            elif mol.get_name() == "si":
                num_si += 1
            elif mol.get_name() in grid:
                num_g += 1

        # Calculate adjustment
        avg = charge/(num_o/2+num_si+num_g)
        q_o = round(self._q_o - avg/2, round_o)
        q_si = round(self._q_si-avg,  round_o)

        # Distribute charge
        for mol in mols:
            if mol.get_name() == "om":
                mol.set_charge(q_o)
            elif mol.get_name() == "si":
                mol.set_charge(q_si)
            elif mol.get_name() in grid:
                mol.set_charge(round(mol.get_charge()-self._q_si+q_si, round_o))

        # Calculate rounding error
        charge_err = sum([mol.get_charge() for mol in mols])
        num_oxy = 100 if num_o > 100 else num_o
        avg_err = charge_err/num_oxy
        q_ox = round(q_o-avg_err, round_o)

        # Randomly distribute balancing oxygens
        if is_rand:
            # Get list of oxygenes
            oxy = {}
            for i in range(len(mols)):
                if mols[i].get_name() == "om":
                    oxy[i] = mols[i]

            # Get random oxygenes and set name and charge
            oxy_r = random.sample([x for x in oxy], num_oxy)
            for o in oxy_r:
                oxy[o].set_name("ox")
                oxy[o].set_short("OX")
                oxy[o].set_charge(q_ox)

            # Push oxy after om
            min_o = min([x for x in oxy])
            max_o = max([x for x in oxy])
            write = []
            write_ox = []
            for i in range(min_o):
                write.append(mols[i])
            for i in range(min_o, max_o+1):
                if i not in oxy_r:
                    write.append(mols[i])
                else:
                    write_ox.append(mols[i])
            write += write_ox
            for i in range(max_o+1, len(mols)):
                write.append(mols[i])

            self._write = write

        # Make global
        self._q = {"om": q_o, "ox": q_ox, "si": q_si}

        self.set_charge(sum([mol.get_charge() for mol in self._write]))


    ############################
    # Public Methods - Analyze #
    ############################
    def _allocation(self, site=None, is_mol=True):
        """Calculate the surface allocation. This is done by first calculating
        the pores surface

        .. math::

            A_{pore} = 2\\pi r\\cdot l_{drill},

        with opening radius :math:`r` and length of the drilling axis
        :math:`l_{drill}`, and the one on the drill sides

        .. math::

            A_{drill} = 2\\cdot\\left(A_{side}-\\pi r^2\\right)

        with block surface on the drilling side :math:`A_{drill}`. Using this
        surface, the allocation rates can be determined by counting the number
        of used sites and free sites from the input binding site list.

        This was done to determine the maximal possible allocation
        :math:`c_{tot}`, the rate for the molecule modification
        :math:`c_{mod}` and the rate for the silanol modification
        :math:`c_{oh}`. For better overview, the relative allocation for the
        molecule modification :math:`c_{rel}` has been calculated

        .. math::

            \\begin{array}{cccc}
            c_{tot}^{pore} = \\dfrac{s_{tot}^{pore}}{A_{pore}},&
            c_{mod}^{pore} = \\dfrac{s_{mod}^{pore}}{A_{pore}},&
            c_{oh}^{pore}  = \\dfrac{s_{oh} ^{pore}}{A_{pore}},&
            c_{rel}^{pore} = \\dfrac{s_{mod}^{pore}}{s_{tot}^{pore}}
            \\end{array}

        with the total number of binding sites :math:`s_{tot}`,
        number of sites used for the modification :math:`s_{mod}`
        and number of sites used for the silanol modification :math:`s_{oh}`.
        The resulting units are

        .. math::

            [c_{i}]=\\frac{\\text{Number of molecules}}{nm^2}
            =\\frac{1}{N_A}\\frac{mol}{nm^2}
            =\\frac{10}{6.022}\\frac{\mu mol}{m^2}.

        Parameters
        ----------
        site : None, list
            Binding site list, None for the public list
        is_mol : bool
            True to calculate allocation in micro molar

        Returns
        -------
        allocation : dict
            Surface allocation values
        """
        # Initialize
        size = self._size
        diam = self._diam

        # User input
        site = self._site if site == None else site
        unit = 10/6.022 if is_mol else 1

        # Calculate surfaces
        surface = {"pore":  size[2]*2*math.pi*(diam/2),
                   "drill": 2*(size[0]*size[1]-math.pi*(diam/2)**2)}

        # Get total number of sites
        count = {}
        count["pore"] = {"tot": sum([1 for x in site if x[4] == 0]),
                         "mod":  sum([1 for x in site if x[4] == 0 and x[3] == 1])}
        count["drill"] = {"tot": sum([1 for x in site if x[4] == 1]),
                          "mod":  sum([1 for x in site if x[4] == 1 and x[3] == 1])}

        # Total allocation
        alloc = {}
        alloc["pore"] = {"tot": count["pore"]["tot"]/surface["pore"]*unit,
                         "mod":  count["pore"]["mod"]/surface["pore"]*unit,
                         "oh": (count["pore"]["tot"]-count["pore"]["mod"])/surface["pore"]*unit,
                         "rel": count["pore"]["mod"]/count["pore"]["tot"]}
        alloc["drill"] = {"tot": count["drill"]["tot"]/surface["drill"]*unit,
                          "mod":  count["drill"]["mod"]/surface["drill"]*unit,
                          "oh": (count["drill"]["tot"]-count["drill"]["mod"])/surface["drill"]*unit,
                          "rel": count["drill"]["mod"]/count["drill"]["tot"]}

        return alloc

    def _rough(self):
        """Calculate the pore roughness. This is normally done by calculating
        the standard deviation of the peaks and valles of a surface.

        In the case of a pore one can visualize pulling the pore appart creating
        a flat surface of the pores interior. Of this surface the roughness is
        determined by calculating the standard deviation of the binding site
        silica peaks and valleys.

        The only difference in the pore is therefore the definition of the axis,
        which is going to be the centeral axis. the mean value :math:`\\bar r`
        of the silica distances :math:`r_i` of silica :math:`i` towards the
        pore center, is calculated by

        .. math::

            \\bar r=\\frac1n\\sum_{i=1}^nr_i

        with the number of silica atoms :math:`n`. This mean value goes into
        the sqareroot roughness calculation

        .. math::

            R_q = \\sqrt{\\frac1n\\sum_{i=1}^n\\|r_i-\\bar r\\|^2}.

        Returns
        -------
        roughness : float
            Pore roughness in nm
        """
        # Initialize
        site = self._site
        center = self._center

        # Run through silica atoms and calculate distances
        r_list = []
        for x in site:
            if x[4] == 0:
                surf = x[1]
                cent = [center[0], center[1], self._data[2][surf]]
                r_list.append(self._length(self._vector(self.pos(surf), cent)))

        # Calculate mean
        r_bar = sum(r_list)/len(r_list)

        # Calculate roughness
        return math.sqrt(sum([(ri-r_bar)**2 for ri in r_list])/len(r_list))

    def _diameter(self):
        """Calculate the resulting pore diameter. This is done by calculating
        the mean value :math:`\\bar r` of the silica distances :math:`r_i` of
        binding site silica :math:`i` towards the pore center

        .. math::

            d=2\\bar r=\\frac2n\\sum_{i=1}^nr_i

        with the number of silica atoms :math:`n`.

        Returns
        -------
        diameter : float
            Pore roughness in nm
        """
        # Initialize
        site = self._site
        center = self._center

        # Run through silica atoms and calculate distances
        r_list = []
        for x in site:
            if x[4] == 0:
                surf = x[1]
                cent = [center[0], center[1], self._data[2][surf]]
                r_list.append(self._length(self._vector(self.pos(surf), cent)))

        # Calculate mean
        r_bar = sum(r_list)/len(r_list)

        # Calculate roughness
        return 2*r_bar


    #########################
    # Public Methods - Edit #
    #########################
    def special(self, mol, sio, num=2, symmetry="point"):
        """Public method of :func:`_special` for a adding specific number of a
        molecule to the inside of the pore. Moreover it is possible to define
        symmetry of the molecules to each other.

        Parameters
        ----------
        mol : Molecule
            Molecule to be added
        sio : list
            List containing the silica id and directional vector of the molecule
        num : int
            Number of molecules to be added
        symmetry : str
            Molecule symmetry orientation to each other

        Examples
        --------
        >>> self.special(Catalysator(),[35,39,37],2)
        >>> self.special(Catalysator(),[35,39,37],2,symmetry="point")
        >>> self.special(Catalysator(),[35,39,37],3,symmetry="random")
        """
        self._special(mol, sio, num, symmetry)

    def couple(self, mol, sio, rate, intent="inside", inp="percent"):
        """Public method of :func:`_couple` for a adding specified rates of
        molecules to the pores inside and outside. The function may contain
        a list for a simultanious coupling inside and outside of the pore or
        just a singular input for a specific placement. In this case the intent
        argument must be provided.

        Parameters
        ----------
        mol : list, mol
            List of two molecules for the inside and outside, entries can be None
        sio : list
            List containing the silica id and directional vector of the molecules
        rate : list, float
            Rate of molecules to be added
        intent : str
            Set to *inside* to add molecules inside the pore and *outside* for
            the outside surface
        inp : str
            Rate type

        Examples
        --------
        >>> self.couple([Epoxi(),Silyl()],[[0,1,2],[0,1,2]],[50,80])
        >>> self.couple([Epoxi2(),None],[[0,1,2],[0,1,2]],[50,80])
        >>> self.couple([None,Silyl()],[[0,1,2],[0,1,2]],[50,80])
        >>> self.couple(Epoxi(),[0,1,2],50)
        >>> self.couple(Silyl(),[0,1,2],80,intent="outside")
        """
        # Process user input
        isList = isinstance(mol, list) and isinstance(sio, list) and isinstance(rate, list)
        if intent == "inside":
            intent = 0
        elif intent == "outside":
            intent = 1
        else:
            print("Wrong intent in coupling function...")
            return

        # Run coupling
        t = utils.tic()

        if isList:
            if mol[0] is not None:
                self._couple(mol[0], sio[0], 0, rate[0], inp=inp)
            if mol[1] is not None:
                self._couple(mol[1], sio[1], 1, rate[1], inp=inp)
        else:
            self._couple(mol, sio, intent, rate, inp=inp)

        self._t_tot += utils.toc(t, "Couple ", self._is_time)

    def coupleProxi(self, mol, sio, orient, rate, inp="percent"):
        """Public method of :func:`_couple_proxi` for a adding specified rates of
        molecules requiring two binding sites to the pore inside.

        Parameters
        ----------
        mol : Molecule
            Molecule to be added
        sio : list
            List containing two silica-oxygen vectors
        orient : list
            Molecule direction vector
        rate : float
            Rate of how many binding sites
        inp : str
            Input type of the rate

        Examples
        --------
        >>> self.coupleProxi(Epoxi(1),[[0,2],[1,3]],[4,6],80)
        >>> self.coupleProxi(Epoxi(1),[[0,2],[1,3]],[4,6],80,inp="precent")
        >>> self.coupleProxi(Epoxi(1),[[0,2],[1,3]],[4,6],10,inp="num")
        """
        t = utils.tic()
        if mol is not None:
            self._couple_proxi(mol, sio, orient, rate, inp=inp)
        self._t_tot += utils.toc(t, "Couple ", self._is_time)

    def finish(self, slx=0, is_rand=True, is_mol=False, is_props=True):
        """Finalize pore by adding siloxan bridges, filling empty bond with
        silanol groups, connecting geminal molecules, removing marked atoms,
        converting the grid to individual molecules and finally moving the pore
        into position. Furthermore the allocation is printed if needed.

        Parameters
        ----------
        slx : float
            Rate of siloxcan bridges
        is_rand : bool
            True to randomly distribute rounding error charge on block oxygen atoms
        is_mol : bool
            True to calculate allocation in micro molar
        is_props : bool
            True to calculate pore properties
        """
        # Start timer
        t = utils.tic()

        # Analysis
        alloc = self._allocation(is_mol=is_mol)  # Calculate allocation
        rough = self._rough()                 # Calculate pore roughness

        # Finalize
        self._siloxan(slx)       # Create silixan bridges
        self._silanol_parallel()  # Fill empty sites with silanol
        self._connect()          # Connect geminal molecules into one
        self._bonding.delete()   # Delete marked atoms
        self._objectify()        # Move all silica and oxygenes into unique mols
        self._sort()             # Sort molecules
        self._excess(is_rand)    # Distribute excess charge
        self._position()         # Position the pore considering pbc
        self._overlap()          # Check for silica or oxygen atoms overlapping

        self._t_tot += utils.toc(t, "Finish ", self._is_time)

        if self._is_time:
            print("----------------------------")
            print("Total   - runtime = "+"%6.3f" % self._t_tot+" s")
            print()

        if is_props:
            # Initialize
            diam = self._diam
            size = self._size
            site = self._site
            drill = self._drill

            # Number of groups
            num_psl = sum([1 for x in site if x[4] == 0 and x[5] is None])
            num_osl = sum([1 for x in site if x[4] == 1 and x[5] is None])
            num_psl_g = sum([1 for x in site if x[4] == 0 and x[5] is not None])
            num_osl_g = sum([1 for x in site if x[4] == 1 and x[5] is not None])

            # Volume and surface
            vol = math.pi*diam**2/4*size[2]
            surf_p = math.pi*diam*size[2]
            surf_o = size[0]*size[1]-math.pi*diam

            # Hydroxylation
            hydro_p = (num_psl+num_psl_g*2)/surf_p
            hydro_o = (num_osl+num_osl_g*2)/surf_o/2

            # Dimensions
            if drill == "x":
                size_str = "%4.1f" % size[2]+" "+"%4.1f" % size[1]+" "+"%4.1f" % size[0]
            elif drill == "y":
                size_str = "%4.1f" % size[0]+" "+"%4.1f" % size[2]+" "+"%4.1f" % size[1]
            elif drill == "z":
                size_str = "%4.1f" % size[0]+" "+"%4.1f" % size[1]+" "+"%4.1f" % size[2]

            # Create dictionary
            props = {
                "D":   ["["+size_str+"]", "nm"],
                "d":   ["%4.2f" % diam,    "nm"],
                "Rq":  ["%5.3f" % rough,   "nm"],
                "S":   ["%6.2f" % surf_p,  "nm^2"],
                "V":   ["%6.2f" % vol,     "nm^3"],
                "hi":  ["%4.2f" % hydro_p, "OH/nm^2"],
                "ho":  ["%4.2f" % hydro_o, "OH/nm^2"],
                "hio": ["%5.3f" % (hydro_p/hydro_o), ""],
                "Ni":  ["%4i" % num_psl,  "#"],
                "Nig": ["%4i" % num_psl_g, "#"],
                "Nii": ["%5.2f" % (num_psl/num_psl_g), ""],
                "No":  ["%4i" % num_osl, " #"],
                "Nog": ["%4i" % num_osl_g, "#"],
                "Noo": ["%5.2f" % (num_osl/num_osl_g), ""],
                "ai":  ["%5.2f" % (alloc["pore"]["rel"]*100), "%"],
                "ao":  ["%5.2f" % (alloc["drill"]["rel"]*100), "%"],
                "aio": ["%5.3f" % (alloc["pore"]["rel"]/alloc["drill"]["rel"]), ""] if alloc["drill"]["rel"] > 0 else ["%5.3f" % 0, ""],
                "C":   ["%7.5e" % self.get_charge(), "C"],
                "t":   ["%5.2f" % self._t_tot,      "s"]
            }
            self._props = props

            # Print allocation
            unit = "[mumol/m^2]" if is_mol else "[#/nm^2]"
            print("Surface allocation in "+unit)
            print("Unmodified Pore  - "+"%6.3f" %
                  alloc["pore"]["tot"] + ",  Out - "+"%6.3f" % alloc["drill"]["tot"])
            print("Modified   Pore  - "+"%6.3f" %
                  alloc["pore"]["mod"] + ",  Out - "+"%6.3f" % alloc["drill"]["mod"])
            print("Silanol    Pore  - "+"%6.3f" %
                  alloc["pore"]["oh"] + ",  Out - "+"%6.3f" % alloc["drill"]["oh"])
            print("Relative   Pore  - "+"%6.3f" %
                  (alloc["pore"]["rel"]*100)+"%, Out - "+"%6.3f" % (alloc["drill"]["rel"]*100)+"%")
            print()

            print("Properties")
            print("Dimensions            - "+props["D"][0]+" "+props["D"][1])
            print("Diameter              - "+props["d"][0]+" "+props["d"][1])
            print("Roughness             - "+props["Rq"][0]+" "+props["Rq"][1])
            print("Surface               - "+props["S"][0]+" "+props["S"][1])
            print("Volume                - "+props["V"][0]+" "+props["V"][1])
            print("Hydroxylation in/out  - "+props["hi"][0]+"/" +
                  props["ho"][0]+" "+props["hi"][1]+" = "+props["hio"][0])
            print("Number of pore SL/SLG - "+props["Ni"][0]+"/" +
                  props["Nig"][0]+" "+props["Ni"][1]+" = "+props["Nii"][0])
            print("Number of out  SL/SLG - "+props["No"][0]+"/" +
                  props["Nog"][0]+" "+props["No"][1]+" = "+props["Noo"][0])
            print("Relative allocation   - "+props["ai"][0]+"/" +
                  props["ao"][0]+" "+props["ai"][1]+" = "+props["aio"][0])
            print("Excess charge         - "+props["C"][0]+" "+props["C"][1])
            print("Total time            - "+props["t"][0]+" "+props["t"][1])
            print()


    ############################
    # Public methods - Analyze #
    ############################
    def volume(self):
        """This function calculates the available volume in the pore system.

        Returns
        -------
        volume : float
            System volume
        """
        # Initialize
        size = self._size

        # Calculate volumes
        res = self._res
        for i in range(self._dim):
            if not i == 2:
                res *= size[i]

        pore = math.pi*(self._diam/2)**2

        return 2*res+pore

    ##################
    # Setter Methods #
    ##################
    def set_grid(self, grid):
        """Set the grid molecule list.

        Parameters
        ----------
        grid : list
            List of grid molecule identifier
        """
        self._grid = grid if isinstance(grid, list) else [grid]


    ##################
    # Getter Methods #
    ##################
    def get_time(self):
        """Return the total creation time.

        Returns
        -------
        time : float
            Used pore creation time
        """
        return self._t_tot

    def get_size(self):
        """Return the pore size.

        Returns
        -------
        size : list
            Pore dimension
        """
        return self._size

    def get_diam(self):
        """Return the pore diameter.

        Returns
        -------
        diam : float
            Pore diameter
        """
        return self._diam

    def get_reservoir(self):
        """Return the reservoir length.

        Returns
        -------
        res : float
            Reservoir length
        """
        return self._res

    def get_grid(self):
        """Return the grid molecules.

        Returns
        -------
        grid : list
            List of grid molecule identifier
        """
        return self._grid

    def get_q(self):
        """Return the adjusted grid charges.

        Returns
        -------
        q : dict
            Dictionary of grid molecules and their charges
        """
        return self._q

    def get_drill(self):
        """Return the drilling axis.

        Returns
        -------
        drill : str
            Drilling axis
        """
        return self._drill

    def get_gap(self):
        """Return the gap vector.

        Returns
        -------
        gap : list
            Gap vector
        """
        return self._gap

    def get_center(self):
        """Return the pore focal point.

        Returns
        -------
        center : list
            Pore focal point
        """
        return self._center

    def get_props(self):
        """Return the final pore properties.

        Returns
        -------
        props : dict
            Dictionary containing properties of the generated pore
        """
        return self._props
