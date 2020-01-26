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
    silicon pores. Note that this class can therefore use all methods of the
    Molecule class.

    At the beginning a silicon-oxygen grid is generated from the smallest
    possible grid unit by duplication. Using the verlet list
    :class:`porems.verlet.Verlet` and bonding classes
    :class:`porems.bonding.Bonding`, this silicon grid is prepared rapidly
    and binding sites are provided.

    Here methods are provided for adding molecule groups to the
    surface of the generated pore. Furthermore, a method is added to fill all
    unused binding sites with silanol and siloxane groups.

    All bound molecules are added to a global molecule list, in order to
    preserve individual molecule properties.

    Parameters
    ----------
    size : list, float
        Size of the silicon-oxygen-grid, can be a float, if same size in all
        dimensions
    diam : float
        Pore diameter
    drill : string
        Axis for the pore drilling
    res : float, optional
        Value to extend and translate the pore by on the drill side creating
        reservoirs
    vs : float, optional
        Verlet-list box size
    proxi_dist : float
        Minimal distance between two groups attached to the surface
    is_pbc : bool, optional
        True if periodic boundary conditions are needed
    is_time : bool, optional
        True to print the used time for the single steps

    Examples
    --------
    Following example generates a pore with TMS molecules on the inside and
    outside surface

    .. code-block:: python

        from porems.pore import Pore
        from porems.essentials import TMS

        pore = Pore([10, 10, 10], 6, "z", 5.5)

        pore.attach(TMS(), [0, 1], [1, 2], 0, 4.3, inp="molar")
        pore.attach(TMS(), [0, 1], [1, 2], 1, 60, inp="percent")

        pore.finalize()
    """
    def __init__(self, size, diam, drill, res=1, vs=0.4, proxi_dist=0.5, is_pbc=True, is_time=False):
        # Call super class
        super(Pore, self).__init__()

        # Initialize
        self._name = "pore"
        self._size = size if isinstance(size, list) else [size for s in range(self._dim)]
        self._diam = diam-{"x": 0.55, "y": 0.5, "z": 0.5}[drill]
        self._is_pbc = is_pbc
        self._is_time = is_time
        self._res = res
        self._vs = vs
        self._o_proxi = proxi_dist
        self._drill = drill

        # Define bond lengths
        self._si_grid = 0.155
        self._si_oh = 0.164
        self._oh = 0.098

        # Define silicon block information
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
        self._t_tot = {}

        # Properties
        self._props = {}
        self._is_props = False
        self._slx = 0

        # Create silicon pore
        t = utils.tic()
        self._block(0, self._foundation())             # Build block
        self._orientation()                            # Rotate drill axis
        self.translate(self._gap)                      # Translate gap
        self._center = self._focal()                   # Find focal point
        self._t_tot["Build"] = utils.toc(t, "Build   ", is_time)

        # Verlet and bonding
        t = utils.tic()
        self._verlet = Verlet(self, vs, is_pbc)        # Create verlet boxes
        self._bonding = Bonding(self._verlet)          # Create bond matrix
        self._t_tot["Matrix"] = utils.toc(t, "Matrix  ", is_time)

        # Prepare pore
        t = utils.tic()
        self._box = self.get_box()                     # Set box
        self._bonding.attach()                         # Prepare sites
        self._bonding.drill(self._center, self._diam)  # Drill pore
        self._bonding.prepare()                        # Prepare pore surface
        self._t_tot["Prepare"] = utils.toc(t, "Prepare ", is_time)

        # Find binding sites
        t = utils.tic()
        self.zero()
        self._box = self.get_box()                     # Reset box size
        self._bind()                                   # Create site array
        self._proxi()                                  # Find sites in proximity
        self._t_tot["Binding"] = utils.toc(t, "Binding ", is_time)

        # Calculate properties
        t = utils.tic()
        self._calc_props()                             # Recalculate properties
        self._diam = self._props["Diameter"]           # Set new diameter
        self._t_tot["Props"] = utils.toc(t, "Props   ", is_time)


    ############
    # Building #
    ############
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
        block.delete([50, 48, 36, 40, 41]) # y-Axis
        block.delete([45, 20, 19, 21, 22, 23, 24, 25, 17, 18]) # z-Axis

        # Move to zero
        block.zero()

        return block

    def _block(self, dim, block):
        """Recursively duplicate and translate a given molecule block in all
        given dimensions. The duplication stops if the next added block would
        create the pore longer than specified in the constructor.

        Parameters
        ----------
        dim : integer
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


    #################
    # Binding sites #
    #################
    def _bind(self):
        """Create binding site matrix. This matrix :math:`\\boldsymbol{B}` is
        mostly generated by the bonding class :class:`poresms.bonding.Bonding`

        .. math::

            \\boldsymbol{B}=\\begin{bmatrix}
                o_0&s_0&p_0&u_0&t_0&g_0&m_0\\\\
                o_1&s_1&p_1&u_1&t_1&g_1&m_1\\\\
                \\vdots\\\\
                o_b&s_b&p_b&u_b&t_0&g_b&m_b\\\\
            \\end{bmatrix}

        with binding site number :math:`b` and entries

        0. Oxygen atom id :math:`o_{k=0,\\dots,b}`
        1. Silica atom id :math:`s_k`
        2. List of pointers :math:`p_k` to binding sites that are in proximity
        3. State :math:`u_k`: available - 0, used - 1, is in proximity - 2
        4. Type :math:`t_k`: inside the pore - 0, outsite the pore - 1
        5. Pointer :math:`g_k` to geminal binding site
        6. Molecule position :math:`m_k` in the global molecule list
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
        by looking for oxygen atoms that are near each other.

        Fill the result in the binding site matrix in the proximity entry.
        """
        # Initialize
        site = self._site

        # Create local verlet list
        atoms = [o[0] for o in site]
        verlet = Verlet(self._temp(atoms), self._o_proxi, self._is_pbc)

        # Find oxygens in proximity
        box_list = [i for i in range(len(verlet.get_box()[1]))]
        oo = verlet.find_parallel(box_list, ["O", "O"], 0.0, self._o_proxi)
        oo_col = utils.column(oo)

        # Append proxi list to sites
        for i in range(len(site)):
            os = site[i]
            o = oo[oo_col[0].index(i)][1] if i in oo_col[0] else []

            os[2].extend(o)

    def _random(self, site_type, rate, inp, counter, is_proxi=False):
        """Get a random binding site list for a given percentage of the
        available binding sites.

        If a site is available (state is 0), it is added to the output-list and
        the state entry for this site is set to 1.
        Then the states of all binding sites in the proximity list are set to 2.

        In case the random site is not available, this run is repeated.
        A break-counter prevents the occurrence of an endless loop.

        The ``site_type`` determines whether the search is done on the inside of
        the pore **0** or on the outside surface **1**.

        Currently following rate types for input ``inp`` are available

        * **num** - For a specific number of binding sites
        * **percent** - For a percentage of the total number of binding sites on the selected surface
        * **molar** - For a molar value dependent on the surface area of the selected surface

        If ``is_proxi`` is set to True, then binding sites nearest to each other
        are determined. the return value then consists of lists containing the
        ids of both sites. Geminal pairs are not allowed, since creating the
        topology is too complex. This however might be added in the future.

        Parameters
        ----------
        site_type : integer
            Type of the sites inside-0, outside-1
        rate : float
            Rate of how many binding sites
        inp : string
            Input type of the rate
        counter : integer
            Number of attempts before breaking the loop
        is_proxi : bool, optional
            True to search for site pairs in proximity

        Returns
        -------
        ran_list : list
            List of pointers to random binding sites
        """
        # Initialize
        site = copy.deepcopy(self._site)
        site_len = utils.column(site)[4].count(site_type)
        id_list = [i for i in range(len(site))]

        # Calculate random element number
        if inp == "num":
            ran_num = int(rate)
        elif inp == "percent":
            ran_num = site_len*rate
            ran_num /= 100
            ran_num = math.ceil(ran_num)
        elif inp == "molar":
            ran_num = rate*6.022/10*self._props["Surface"][site_type]
            ran_num = math.ceil(ran_num)

        # Get list of all closest partners for dual attachment
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

            # Check binding site for dual attachment
            if is_proxi:
                is_app = (is_proxi and site[rand][4] == site_type and
                          site[rand][3] == 0 and site[site_min[rand][0]][3] == 0 and
                          site[rand][5] == None and site[site_min[rand][0]][5] == None)
            # Check binding site for singular attachment
            else:
                is_app = site[rand][4] == site_type and site[rand][3] == 0

            # Append to return list
            if is_app:
                # Reset counter
                count = 0

                # Append
                ran_list.append([rand, site_min[rand][0]] if is_proxi else rand)

                # Set site and proximity to used
                site = self._close([rand, site_min[rand][0]] if is_proxi else rand, site)
            else:
                # Add counter
                count += 1

        return ran_list

    def _close(self, sites, site_list=None):
        """This function closes the specified site by setting the used value
        to 1 and the used value of all sites in proximity to 2 in the given
        site list.

        Parameters
        ----------
        sites : integer, list
            identifiers of sites to be closed
        site_list : list, None, optional
            List for closing off binding sites, leave None for global list

        Returns
        -------
        site : list
            Processed site list
        """
        # Initialize
        site = self._site if site_list is None else site_list

        # Process input
        sites = [sites] if isinstance(sites, int) else sites

        # Close sites
        for x in sites:
            site[x][3] = 1
            for prox in site[x][2]:
                if not site[prox][3] == 1:
                    site[prox][3] = 2

        return site

    def _add_mol_list(self, mol_list):
        """Add given molecule list to global molecule list and add the position
        in this list to the corresponding binding site. The input has the format

        * [binding_site_idx, molecule]
        * [[binding_site_idx_1, molecule_1], [binding_site_idx_2, molecule_2], ...]

        Parameters
        ----------
        mol_list : list
            Molecules to be added to the global list
        """
        # Process input
        if not isinstance(mol_list[0],list):
            mol_list = [mol_list]

        # Add to molecule list
        for mol in mol_list:
            self._site[mol[1]][6] = len(self._mol_list)
            self._mol_list.append(mol[0])


    #######################
    # Molecule Attachment #
    #######################
    def attach(self, mol, si_o, orient, site_type, rate, inp="percent", sites=None, counter=1000,is_rotate=True):
        """Add a molecule to the pore binding site.

        The input ``si_o`` defines the silicon atom id and ``orient`` consists of
        two atom ids that determine which molecule bond should be oriented
        towards the centre.

        The binding sites are chosen randomly using function :func:`_random`.
        Available input types are noted in the mentioned function.

        First the molecule is rotated around its directional vector randomly to
        simulate a state closer to reality. Second the molecule is rotated so
        that its directional vector and the binding site vector match. One final
        rotation is applied for molecules on the inside, so the directional
        vector points towards the centre of the pore.

        The molecule is then moved in order to place its silicon atom on the
        silicon of the binding site. The oxygen is then placed on top of the
        binding site oxygen, in order to keep the sites geometry.

        Finally, the atoms of the binding sites are added to the list of atoms
        to be deleted at the end of the construction.

        These molecules are then added to the global molecule list.

        Parameters
        ----------
        mol : Molecule
            Molecule to be added
        si_o : list
            Silicon and oxygen id of the grid atoms of the molecule
        orient : list
            Two atom ids that determine which molecule bond should be oriented
            towards the centre
        site_type : integer
            Type of the sites **inside-0**, **outside-1**
        rate : float
            Rate of how many binding sites
        inp : string, optional
            Input type of the rate
        sites : integer, list, optional
            List index of binding sites
        counter : integer, optional
            Number of attempts before breaking the loop
        is_rotate : bool, optional
            True to randomly rotate molecules on the surface

        Returns
        -------
        mol_list : list
            List of Molecule objects and corresponding site ids
        remove_list : list
            List of atom ids to be deleted from the pore

        Examples
        --------
        .. code-block:: python

            pore.attach(TMS(), [0, 1], [1, 2], 0, 50)
            pore.attach(TMS(), [0, 1], [1, 2], 0, 50, inp="percent")
            pore.attach(TMS(), [0, 1], [1, 2], 1, 2.4, inp="molar")
        """
        # Stop time
        t = utils.tic()

        # Process user input
        sites = [sites] if isinstance(sites, int) else sites

        # Get random list or given list
        site_list = self._random(site_type, rate, inp, counter) if sites is None else sites

        # Molecule orientation vector
        vec_m = mol.bond(orient[0], orient[1])[0]

        # Centre vector independent of z-axis
        center = self._center[:-1]

        # Add molecule
        mol_list = []
        remove_list = []
        for i in sorted(site_list, reverse=True):
            # Initialize
            site = self._site[i]

            # Copy molecule
            temp = copy.deepcopy(mol)

            # Randomly rotate molecule around orientational axis
            if is_rotate: temp.rotate(vec_m, random.randint(1, 180))

            # Set vector of silicon to oxygen atom of binding site
            vec_o = self._vector(site[1], site[0])

            # Adjust molecule vector to binding site vector
            temp.rotate(self._cross(vec_m, vec_o), self.angle(vec_m, vec_o))

            # Adjust molecule to face central axis on the inside
            if site[4] == 0:
                # Get center axis position perpendicular to binding site
                pos_c = center[:]+[self.pos(site[0])[self._dim-1]]

                # Set vector of binding site silicon atom to central position
                vec_c = self._vector(self.pos(site[1]), pos_c)

                # Rotate molecule towards centre
                temp.rotate(self._cross(vec_c, vec_o), -self.angle(vec_c, vec_o))

            # Adjust orientation perpendicular to surface on the outside
            elif site[4] == 1:
                # Set vector towards or away from z-axis depending on which side
                vec_z = [x if self.pos(site[0])[2]>self._size[2]/2 else -1*x for x in self._axis("z")]

                # Rotate molecule perpendicular to surface
                temp.rotate(self._cross(vec_z, vec_o), -self.angle(vec_z, vec_o))

            # Move molecule to position of binding site originating from the oxygen atom
            temp.move(si_o[1], self.pos(site[0]))

            # Move molecule silicon atom to binding site position
            temp.put(si_o[0], self.pos(site[1]))

            # Check if geminal and set names
            if site[5] is not None:
                temp.set_name(temp.get_name()+"g")
                temp.set_short(temp.get_short()+"G")

            # Remove obsolete atoms
            remove_list.extend([site[0], site[1]])

            # Add to molecule list
            if sites == None:
                self._add_mol_list([temp, i])

            mol_list.append([temp, i])

            # Close binding sites
            self._close(i)

        # Finalize
        if not mol.get_name()=="sl":
            # Remove obsolete atoms
            self._bonding.remove(remove_list)

            # Initialize allocation
            key_name = mol.get_short()
            location = "in" if site_type==0 else "out"

            # Add to time dictionary
            self._t_tot["Attach_"+key_name+"_"+location] = utils.toc(t, "Attach  ", self._is_time)

            # Add number of bonds to allocation
            if not key_name in self._props["Allocation"]:
                self._props["Allocation"][key_name] = {0: 0, 1: 0}

            self._props["Allocation"][key_name][site_type] += len(site_list)

        return mol_list, remove_list

    def special(self, mol, si_o, orient, num, symmetry):
        """Add a molecule in a specified orientation with a specific amount.
        Currently following symmetry orientations to each other are available

        * **random** - No symmetry, random gaussian placement
        * **point** - Point symmetry
        * **mirror** - Mirror symmetry

        Parameters
        ----------
        mol : Molecule
            Molecule to be added
        si_o : list
            Silicon and oxygen id of the grid atoms of the molecule
        orient : list
            Two atom ids that determine which molecule bond should be oriented
            towards the centre
        num : integer
            Number of molecules to be added
        symmetry : string
            Molecule symmetry orientation to each other

        Examples
        --------
        .. code-block:: python

            pore.special(TMS(), [0, 1], [1, 2], 2)
            pore.special(TMS(), [0, 1], [1, 2], 2, symmetry="point")
            pore.special(TMS(), [0, 1], [1, 2], 3, symmetry="random")
        """
        # Process input
        if symmetry not in ["random", "point", "mirror"]:
            print("Symmetry type not supported...")
            return

        elif not symmetry == "random" and not num == 0:
            # Calculate geometrical positions
            pos_list = []
            length = self._center[2]*2
            dist = length/num
            start = dist/2

            for i in range(num):
                if symmetry == "point":
                    coeff = -1 if i % 2 == 0 else 1
                elif symmetry == "mirror":
                    coeff = 1

                x = self._center[0]+coeff*self._diam/2
                y = self._center[1]
                z = start+dist*i

                pos_list.append([x, y, z])

            # Find nearest site
            site_min = [[0, None] for i in range(num)]
            for i in range(len(self._site)):
                length = [self._length(self._vector(self.pos(self._site[i][0]), pos)) for pos in pos_list]

                for j in range(num):
                    if site_min[j][1] is None or length[j] < site_min[j][1]:
                        site_min[j] = [i, length[j]]

            sites = [x[0] for x in site_min]

        else:
            sites = None

        mol_list, temp = self.attach(mol, si_o, orient, 0, num, inp="num", sites=sites, is_rotate=False)

        self._add_mol_list(mol_list)

    def _silanol(self, sites=None):
        """Convert the remaining binding sites to silanol groups using
        function :func:`attach`.

        Parameters
        ----------
        sites : list, None, optional
            Site list to go through

        Returns
        -------
        attach : dictionary
            List of Molecule objects and ids and List of atom ids to be deleted
            from the pore as a dictionary
        """
        # Initialize
        center = self._center[:-1]

        # Set site range
        site_list = [i for i in range(len(self._site))] if sites is None else sites
        site_list = [x for x in site_list if not self._site[x][3] == 1]

        # Define temporary molecule object
        mol = Molecule()
        mol.add("Si",[0, 0, 0])
        mol.add("O", 0, r=self._si_oh)
        mol.add("H", 1, r=self._oh)
        mol.set_name("sl")
        mol.set_short("SL")
        mol.set_charge(self._q_si+self._q_oh)

        # Run attach method
        mol_list, remove_list = self.attach(mol, [0, 1], [0, 1], 0, len(site_list), inp="num", sites=site_list)

        if sites is None:
            self._add_mol_list(mol_list)
            self._bonding.remove(remove_list)
        else:
            return {"mol_list": mol_list, "remove_list": remove_list}

    def _silanol_parallel(self):
        """Parallelized function :func:`_silanol`.
        """
        # Initialize
        np = mp.cpu_count()

        # Split site list
        site_len = math.floor(len(self._site)/np)
        site_list = []
        for i in range(np):
            if i == np-1:
                site_list.append(list(range(site_len*i, len(self._site))))
            else:
                site_list.append(list(range(site_len*i, site_len*(i+1))))

        # Paralellize
        pool = mp.Pool(processes=np)
        results = pool.map_async(self._silanol, site_list)
        pool.close()

        # Extract results
        mol_list = []
        remove_list = []
        for result in results.get():
            mol_list.extend(result["mol_list"])
            remove_list.extend(result["remove_list"])

        # Add to mol_list and add site - molecule connection
        self._add_mol_list(mol_list)

        # Remove obsolete atoms
        self._bonding.remove(remove_list)

    def attach_dual(self, mol, si_o, orient, rate, inp="percent", counter=1000):
        """Add a molecule to two pore binding sites.

        The ``si_o`` input defines the silicon and oxygen atom pair in two
        separate lists.

        First the molecule is rotated, so that the oxygens are on top of the
        binding site oxygens. Second the molecule is rotated so that the
        directional vector points to the centre of the pore.

        The molecule is then moved, so that its oxygen atom is placed on the
        oxygen of the first binding site then translated, so the molecule is in
        the centre of both sites.

        Finally, the atoms of the binding sites are added to the list of atoms
        to be deleted at the end of the construction.

        Parameters
        ----------
        mol : Molecule
            Molecule to be added
        si_o : list
            List containing two silicon-oxygen vectors
        orient : list
            Molecule direction vector
        rate : float
            Rate of how many binding sites
        inp : string, optional
            Input type of the rate
        counter : integer, optional
            Number of attempts before breaking the loop

        Examples
        --------
        .. code-block:: python

            pore.attach_dual(DualMol(), [[0, 2], [1, 3]], [4, 6], 80)
            pore.attach_dual(DualMol(), [[0, 2], [1, 3]], [4, 6], 80, inp="precent")
            pore.attach_dual(DualMol(), [[0, 2], [1, 3]], [4, 6], 10, inp="num")
        """
        # Get random list
        site_list = self._random(0, rate, inp, counter, is_proxi=True)

        # Molecule and centre vector
        center = self._center[:-1]

        # Stop time
        t = utils.tic()

        # Add molecule
        for i in sorted(site_list, reverse=True):
            # Initialize
            sites = [self._site[i[0]], self._site[i[1]]]

            # Copy molecule
            temp = copy.deepcopy(mol)

            # Rotate towards binding site
            vec_mo = temp.bond(si_o[0][1], si_o[1][1])[0]   # Molecule oxygens
            vec_so = self._vector(sites[0][0], sites[1][0]) # Binding site oxygens
            vec_m = temp.bond(orient[0], orient[1])[0]      # Molecule orientation

            temp.rotate(vec_m, self.angle(vec_mo, vec_so))

            # Rotate towards central axis
            cent_o = [(self.pos(sites[0][0])[x]+self.pos(sites[1][0])[x])/2 for x in range(self._dim)]
            cent_p = center[:]+[cent_o[-1]]

            vec_m = temp.bond(orient[0], orient[1])[0]  # Molecule orientation
            vec_c = self._vector(cent_o, cent_p)        # Binding site center

            temp.rotate(self._cross(vec_m, vec_c), self.angle(vec_m, vec_c))

            # Put oxygens on top of each other
            temp.move(si_o[0][1], self.pos(sites[0][0]))

            # Move to centre
            cent_m = [(temp.pos(si_o[0][1])[x]+temp.pos(si_o[1][1])[x])/2 for x in range(self._dim)]
            temp.translate(self._vector(cent_m, cent_o))

            # Move silicon atoms
            temp.put(si_o[0][0], self.pos(sites[0][1]))
            temp.put(si_o[1][0], self.pos(sites[1][1]))

            # Remove atoms
            self._bonding.remove([sites[0][0], sites[0][1], sites[1][0], sites[1][1]])

            # Add to molecule list
            self._mol_list.append(temp)

            # Close binding sites
            self._close(i)

        # Finalize
        key_name = mol.get_name()+"_0"

        # Add to time dictionary
        self._t_tot["Attach_"+key_name] = utils.toc(t, "Attach  ", self._is_time)

        # Add number of bonds to allocation
        if not key_name in self._props["Allocation"]:
            self._props["Allocation"][key_name] = 0
        self._props["Allocation"][key_name] += len(site_list)*2

    def siloxane(self, rate, inp="percent", counter=1000):
        """Add siloxane bridges to the pore by selecting two binding sites in
        proximity removing both oxygens and placing one at centre of
        the two removed atoms.

        Parameters
        ----------
        rate : float
            Rate of binding sites to be edited
        inp : string, optional
            Rate input type
        counter : integer, optional
            Number of attempts before breaking the loop

        Examples
        --------
        .. code-block:: python

            pore.siloxane(10)
            pore.siloxane(10, inp="precent")
            pore.siloxane(15, inp="num")
            pore.siloxane(2.4, inp="molar")
        """
        # Exit if rate is zero
        if rate == 0:
            return

        # Get random list
        site_list = self._random(0, rate, inp, counter, is_proxi=True)
        self._slx += len(site_list)

        # Molecule and centre vector
        center = self._center[:-1]

        # Define temporary molecule object
        mol = Molecule()
        mol.set_name("slx")
        mol.set_short("SLX")
        mol.add("OM", [0, 0, 0])
        mol.add("O", 0, r=0.1)
        mol.add("O", 0, r=0.1)

        # Add molecule
        for i in sorted(site_list, reverse=True):
            # Initialize
            site = [self._site[i[0]], self._site[i[1]]]

            # Get centre position
            cent_o = [(self.pos(site[0][0])[x]+self.pos(site[1][0])[x])/2 for x in range(self._dim)]
            cent_s = [(self.pos(site[0][1])[x]+self.pos(site[1][1])[x])/2 for x in range(self._dim)]

            # Create molecule
            temp = copy.deepcopy(mol)
            temp.put(0, cent_o)
            temp.translate(self._vector(cent_s, cent_o))
            temp.put(1, self.pos(site[0][0]))
            temp.put(2, self.pos(site[1][0]))

            # Rotate towards central axis
            cent_p = center[:]+[cent_o[-1]]

            vec_m = self._vector(cent_o, temp.pos(0))
            vec_c = self._vector(cent_o, cent_p)

            temp.rotate(self._cross(vec_m, vec_c), self.angle(vec_m, vec_c))

            # Put oxygens on top of each other
            temp.move(1, self.pos(site[0][0]))

            # Remove temporary oxygen atoms
            temp.delete([1, 2])

            # Remove atoms
            self._bonding.remove([site[0][0], site[1][0]])

            # Add to molecule list
            self._mol_list.append(temp)

            # Close binding sites
            self._close(i)


    ###############
    # Final Edits #
    ###############
    def _connect(self):
        """Concatenate two geminal molecules to one molecule.
        """
        # Initialize
        site = self._site
        mol_list = self._mol_list
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
            id_a = 0 if mol_list[g[1]].get_name() == "slg" else 1
            id_b = abs(id_a-1)

            mol_a = mol_list[g[id_a]]
            mol_b = mol_list[g[id_b]]

            num = mol_b.get_num()

            for i in range(num):
                if not mol_b.get_type(i) == "Si":
                    mol_a.add(mol_b.get_type(i), mol_b.pos(i))

            pop.append(g[id_b])

            # Add charge
            mol_list[g[id_a]].set_charge(mol_list[g[id_a]].get_charge()+self._q_oh)

        # Delete from mol_list list
        for d in sorted(pop, reverse=True):
            mol_list.pop(d)

    def _objectify(self, data=None):
        """Move all remaining grid silicon and oxygen atoms to individual
        molecules for writing the structure file and add these molecules to the
        global molecule list.

        Parameters
        ----------
        data_range : list, None, optional
            Datarange to process

        Returns
        -------
        mol_list : list
            Molecule list to be added to global molecule list
        """
        # Define temporary molecule object
        temp_mol = Molecule()
        data_range = [0,len(self._data[0])] if data is None else data

        # Create OM and SI mols
        mol_list = []
        for i in range(data_range[0], data_range[1]):
            temp = copy.deepcopy(temp_mol)

            # Check if oxygen
            is_oxy = self.get_type(i) == "O"

            # Set temporary vars
            atom = "OM" if is_oxy else "SI"
            name = "om" if is_oxy else "si"
            short = "OM" if is_oxy else "SI"
            charge = self._q_o if is_oxy else self._q_si

            # Testing
            if self.get_type(i) == "R":
                atom = "R"
                name = "test"

            # Create unique object
            temp.add(atom, self.pos(i))
            temp.set_name(name)
            temp.set_short(short)
            temp.set_charge(charge)

            # Add to molecule list
            if data is None:
                self._mol_list.append(temp)
            else:
                mol_list.append(temp)

        # Add mols to molecule list
        if data is None:
            self._mol_list.pop(0)
        else:
            return mol_list

    def _objectify_parallel(self):
        """Parallelized function :func:`_objectify`.
        """
        # Initialize
        np = mp.cpu_count()

        # Split site list
        data_len = math.floor(len(self._data[0])/np)
        data_list = []
        for i in range(np):
            if i == np-1:
                data_list.append([data_len*i, len(self._data[0])])
            else:
                data_list.append([data_len*i, data_len*(i+1)])

        # Parallelize
        pool = mp.Pool(processes=np)
        results = pool.map_async(self._objectify, data_list)
        pool.close()

        # Extract results
        mol_list = []
        for result in results.get():
            mol_list.extend(result)

        # Add to molecule list
        for mol in mol_list:
            self._mol_list.append(mol)

        # Remove pore from molecule list
        self._mol_list.pop(0)

    def _position(self):
        """Translate the pore into position and consider the periodic repeat
        gap. Hereby the coordinates are moved by halve the gap distance and box
        is extended by the other halve. Additionally, to prevent molecule
        overlap and periodic movement, the drill side is extended and translated
        by a fixed value, thus, creating reservoirs.
        """
        # Get vector to coordinate zero
        vec = Molecule(inp=[mol for mol in self._mol_list if mol.get_name() in ["si", "om"]]).zero()

        # Position pore
        for mol in self._mol_list:
            mol.translate([vec[i]+self._gap[i]/2 for i in range(self._dim)])
            mol.translate([self._res if i == 2 else 0 for i in range(self._dim)])

        # Extend box
        coord = Molecule(inp=[mol for mol in self._mol_list if mol.get_name() in ["si", "om", "ox"]]).get_box()

        self._box = [coord[i] + self._gap[i]/2 for i in range(self._dim)]

        self._box[2] += self._res

    def _overlap(self):
        """Method for checking of silicon or oxygen atoms are overlapping
        using verlet lists. Duplicate atoms will be printed.
        """
        # Initialize
        temp = Molecule(inp=self._mol_list)

        # Create verlet
        verlet = Verlet(temp, self._vs, True)

        # Check silicon and oxygen
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
        """Adjust grid atoms charges to compensate for excess charges on the
        pore.

        Parameters
        ----------
        is_rand : bool, optional
            True to randomly distribute rounding error charge on block oxygen
            atoms
        """
        # Initialize
        mols = self._mol_list

        # Get excess charge
        charge = 0.0
        num_o = 0
        num_si = 0
        num_g = 0
        round_o = 6
        for mol in mols:
            charge += mol.get_charge()

            if mol.get_name() == "om":
                num_o += 1
            elif mol.get_name() == "si":
                num_si += 1
            else:
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
            else:
                mol.set_charge(round(mol.get_charge()-self._q_si+q_si, round_o))

        # Calculate rounding error
        charge_err = sum([mol.get_charge() for mol in mols])
        num_oxy = 100 if num_o > 100 else num_o
        avg_err = charge_err/num_oxy
        q_ox = round(q_o-avg_err, round_o)

        # Randomly distribute balancing oxygens
        if is_rand:
            # Get list of oxygens
            oxy = {}
            for i in range(len(mols)):
                if mols[i].get_name() == "om":
                    oxy[i] = mols[i]

            # Get random oxygens and set name and charge
            oxy_r = random.sample([x for x in oxy], num_oxy)
            for o in oxy_r:
                oxy[o].set_name("ox")
                oxy[o].set_short("OX")
                oxy[o].set_charge(q_ox)

            # Push oxy after om
            min_o = min([x for x in oxy])
            max_o = max([x for x in oxy])
            mol_list = []
            mol_list_ox = []
            for i in range(min_o):
                mol_list.append(mols[i])
            for i in range(min_o, max_o+1):
                if i not in oxy_r:
                    mol_list.append(mols[i])
                else:
                    mol_list_ox.append(mols[i])
            mol_list += mol_list_ox
            for i in range(max_o+1, len(mols)):
                mol_list.append(mols[i])

            self._mol_list = mol_list

        # Make global
        self._q = {"om": q_o, "ox": q_ox, "si": q_si}

        self.set_charge(sum([mol.get_charge() for mol in self._mol_list]))

    def _sort(self):
        """Sort molecules in order to prevent multiple molecule appearances in
        the structure file.
        """
        # Add molecules to a dictionary
        mol_dict = {}
        for mol in self._mol_list:
            if mol.get_name() not in mol_dict:
                mol_dict[mol.get_name()] = []

            mol_dict[mol.get_name()].append(mol)

        # Create new molecule list
        mol_list = []
        priority = ["si","om","ox","slx","sl","slg"]
        for mol_name in priority:
            if mol_name in mol_dict:
                mol_list += mol_dict[mol_name]

        for mol_name in mol_dict:
            if mol_name not in priority:
                mol_list += mol_dict[mol_name]

        # Replace global list
        self._mol_list = mol_list


    ###########
    # Analyse #
    ###########
    def _calc_props(self):
        """This function calculates all necessary properties of the system and
        returns a dictionary. Note that most of the properties can only be
        calculated after the pore is finalized in function :func:`finalize`.

        The output dictionary contains following properties which will be
        explained further down

        * **Allocation** - Dictionary of all molecule allocations in molar and percent
        * **Charge** - Excess charge
        * **Diameter** - Generated pore diameter
        * **Dimension** - Generated Dimensions
        * **Generation_Time** - Dictionary of all time steps
        * **Hydroxylation** - Hydroxylation dictionary inside and outside of the pore
        * **Roughness** - Surface roughness
        * **Silanol_Geminal** - Dictionary containing the number of silanol and geminal silanol molecules inside and outside the pore
        * **Surface** - Dictionary of inner and outer pore surface
        * **System_Size** - Total system size including reservoirs
        * **Volume** - Inner pore volume

        The **roughness** is calculated as the standard deviation of the peaks
        and valleys of a surface.

        In the case of a pore one can visualize pulling it apart, creating a
        flat surface out of the interior. The roughness is thus determined by
        calculating the standard deviation of the binding site silicon peaks and
        valleys.

        The only difference in the pore is therefore the definition of the axis,
        which is going to be the central axis. The mean value :math:`\\bar r`
        of the silicon distances :math:`r_i` of silicon :math:`i` towards the
        pore centre, is calculated by

        .. math::

            \\bar r=\\frac1n\\sum_{i=1}^nr_i

        with the number of silicon atoms :math:`n`. This mean value is used in
        the square root roughness calculation

        .. math::

            R_q = \\sqrt{\\frac1n\\sum_{i=1}^n\\|r_i-\\bar r\\|^2}.


        The **diameter** is the mean value :math:`\\bar r` of the silicon
        distances :math:`r_i` of binding site silicon :math:`i` towards the pore
        centre

        .. math::

            d=2\\bar r=\\frac2n\\sum_{i=1}^nr_i.

        The **surface** inside the pore is calculated by

        .. math::

            A_\\text{inside}=2\\pi r\\cdot l_\\text{drill},

        with opening radius :math:`r` and length of the drilling axis
        :math:`l_\\text{drill}`, the total outside surface of the pore system by

        .. math::

            A_\\text{outside}=2\\cdot\\left(A_\\text{side}-\\pi r^2\\right)

        with block surface on the drilling side :math:`A_\\text{side}`,
        and the **volume** by

        .. math::

            V_\\text{pore}=\\pi r^2\\cdot l_\\text{drill}.

        Using these surfaces, the **allocation** rates can be determined by
        counting the number of used sites and free sites from the binding
        site list. Thus, it is possible to determine the maximal possible
        allocation :math:`c_\\text{tot}` also called hydroxylation, the rate for
        the molecule modification :math:`c_\\text{mod}` and the rate for the
        silanol modification :math:`c_\\text{oh}`. For better overview, the
        relative allocation for the molecule modification :math:`c_\\text{rel}`
        has been calculated

        .. math::

            \\begin{array}{cccc}
            c_\\text{tot} = \\dfrac{s_\\text{tot}}{A},&
            c_\\text{mod} = \\dfrac{s_\\text{mod}}{A},&
            c_\\text{oh}  = \\dfrac{s_\\text{oh} }{A},&
            c_\\text{rel} = \\dfrac{s_\\text{mod}}{s_\\text{tot}}
            \\end{array}

        with the corresponding surface :math:`A` for either inside or outside
        the pore, the total number of binding sites :math:`s_\\text{tot}`,
        number of sites used for the modification :math:`s_\\text{mod}`
        and number of sites used for the silanol modification
        :math:`s_\\text{oh}`. The resulting units are

        .. math::

            [c_{i}]=\\frac{\\text{Number of molecules}}{nm^2}
            =\\frac{1}{N_A}\\frac{mol}{nm^2}
            =\\frac{10}{6.022}\\frac{\\mu mol}{m^2}.
        """
        # Initialize
        props = self._props

        # Calculate generation properties
        if not self._is_props:
            # Run through binding site silicon atoms and calculate distances
            r_list = []
            for site in self._site:
                if site[4] == 0:
                    si_id = site[1]
                    center = [self._center[0], self._center[1], self._data[2][si_id]]
                    r_list.append(self._length(self._vector(self.pos(si_id), center)))

            # Calculate mean
            r_bar = sum(r_list)/len(r_list) if len(r_list)>0 else 0

            # Set dimension order
            dim_order = {"x": [2, 1, 0], "y": [0, 2, 1], "z": [0, 1, 2]}

            # Calculate properties
            props["Roughness"] = math.sqrt(sum([(ri-r_bar)**2 for ri in r_list])/len(r_list)) if r_bar>0 else 0
            props["Diameter"] = 2*r_bar if r_bar>0 else 0
            props["Surface"] = {0: self._size[2]*2*math.pi*(props["Diameter"]/2),
                                1: 2*(self._size[0]*self._size[1]-math.pi*(props["Diameter"]/2)**2)}
            props["Volume"] = self._size[2]*math.pi*(props["Diameter"]/2)**2
            props["Dimension"] = [self._size[dim_order[self._drill][i]] for i in range(3)]
            props["Allocation"] = {}

        # Calculate finished properties
        else:
            # Molar conversion
            molar = 10/6.022

            # Count number of groups
            num_i = sum([1 for x in self._site if x[4] == 0 and x[5] is None])
            num_o = sum([1 for x in self._site if x[4] == 1 and x[5] is None])
            num_i_g = sum([1 for x in self._site if x[4] == 0 and x[5] is not None])/2
            num_o_g = sum([1 for x in self._site if x[4] == 1 and x[5] is not None])/2

            num_tot = [num_i+num_i_g-2*self._slx, num_o+num_o_g]

            # Calculate hydroxylation
            hydro = [num_tot[0]/props["Surface"][0], num_tot[1]/props["Surface"][1]]

            # Calculate allocation
            num_oh = [num_i+num_i_g, num_o+num_o_g]

            for mol in props["Allocation"]:
                mol_data = {}
                for site_type in props["Allocation"][mol]:
                    location = "in" if site_type==0 else "out"

                    mod = props["Allocation"][mol][site_type]/props["Surface"][site_type]
                    rel = props["Allocation"][mol][site_type]/num_tot[site_type]

                    mol_data[location] = [props["Allocation"][mol][site_type], mod, mod*molar]
                    # mol_data[location] = [rel, rel*molar] # Relative allocation

                    # Count oh
                    num_oh[site_type] -= props["Allocation"][mol][site_type]

                props["Allocation"][mol] = mol_data

            # OH allocation
            oh = [num_oh[0]/props["Surface"][0], num_oh[1]/props["Surface"][1]]
            props["Allocation"]["OH"] = {"in": [int(num_oh[0]), oh[0], oh[0]*molar], "out": [int(num_oh[1]), oh[1], oh[1]*molar]}

            # Siloxane allocation
            slx = self._slx/props["Surface"][0]
            props["Allocation"]["SLX"] = {"in": [self._slx, slx, slx*molar], "out": [0, 0.0, 0.0]}

            # Calculate properties
            props["Charge"] = self._charge
            props["System_Size"] = self._box
            props["Generation_Time"] = self._t_tot
            props["Hydroxylation"] = {"in": [hydro[0],hydro[0]*molar], "out": [hydro[1],hydro[1]*molar]}
            props["Silanol_Geminal"] = {"in": num_i, "out": num_o, "in_g": num_i_g, "out_g": num_o_g}


    ################
    # Finalization #
    ################
    def finalize(self, is_rand=True):
        """Finalize pore by filling empty bond with
        silanol groups, connecting geminal molecules, removing marked atoms,
        converting the grid to individual molecules and finally moving the pore
        into position.

        Parameters
        ----------
        is_rand : bool, optional
            True to randomly distribute rounding error charge on block oxygen
            atoms
        """
        # Fill empty sites with silanol
        t = utils.tic()
        self._silanol_parallel()
        self._t_tot["Silanol"] = utils.toc(t, "Silanol ", self._is_time)

        # Connect geminal molecules into one
        t = utils.tic()
        self._connect()
        self._t_tot["Connect"] = utils.toc(t, "Connect ", self._is_time)

        # Delete marked atoms
        t = utils.tic()
        self._bonding.delete()
        self._t_tot["Delete"] = utils.toc(t, "Delete  ", self._is_time)

        # Move all silicon and oxygens into unique mols
        t = utils.tic()
        self._objectify_parallel()
        self._t_tot["Objects"] = utils.toc(t, "Objects ", self._is_time)

        # Distribute excess charge
        t = utils.tic()
        self._excess(is_rand)
        self._t_tot["Excess"] = utils.toc(t, "Excess  ", self._is_time)

        # Sort molecules
        t = utils.tic()
        self._sort()
        self._t_tot["Sort"] = utils.toc(t, "Sort    ", self._is_time)

        # Position the pore considering pbc
        t = utils.tic()
        self._position()
        self._t_tot["Position"] = utils.toc(t, "Position", self._is_time)

        # Check for silicon or oxygen atoms overlapping
        t = utils.tic()
        self._overlap()
        self._t_tot["Overlap"] = utils.toc(t, "Overlap ", self._is_time)

        # Allow properties calculation
        self._is_props = True
        self._calc_props()

        if self._is_time:
            print("----------------------------")
            print("Total   - runtime = "+"%6.3f" % sum([self._t_tot[x] for x in self._t_tot])+" s")
            print()

        # Print finalization note
        print("Pore is finalized. Attachment functions should not be called anymore.")

    def table_props(self, decimals=3):
        """This functions creates readable pandas DataFrames of all properties
        data. Available tables are

        * **props** - Pore properties
        * **alloc** - Surface allocation
        * **time** - Time consumption

        Parameters
        ----------
        decimals : integer
            Number of decimals to be rounded to

        Returns
        -------
        tables : dictionary
            Dictionary of pandas DataFrames
        """
        # Import
        import pandas as pd

        # Initialize
        tables = {}
        props = self._props
        form = "%."+str(decimals)+"f"

        # Properties table
        data_props = {}
        data_props["Dimension"] = "["+form%props["Dimension"][0]+", "+form%props["Dimension"][1]+", "+form%props["Dimension"][2]+"]"
        data_props["Diameter"] = form%props["Diameter"]
        data_props["Roughness"] = form%props["Roughness"]
        data_props["Surface"] = form%props["Surface"][0]
        data_props["Volume"] = form%props["Volume"]
        data_props["Hydroxylation"] = form%props["Hydroxylation"]["in"][0]+"/"+form%props["Hydroxylation"]["out"][0]
        data_props["Silanol_In"] = "%i"%props["Silanol_Geminal"]["in"]+"/"+"%i"%props["Silanol_Geminal"]["in_g"]
        data_props["Silanol_Out"] = "%i"%props["Silanol_Geminal"]["out"]+"/"+"%i"%props["Silanol_Geminal"]["out_g"]
        data_props["Time"] = form%sum([props["Generation_Time"][x] for x in props["Generation_Time"]])

        tables["props"] = pd.DataFrame.from_dict(data_props, orient="index", columns={self._drill+"-Axis - "+"%.0f"%self._diam+"nm"})

        # Allocation table
        data_alloc = {}
        for mol in props["Allocation"]:
            data_alloc[mol] = {}
            for location in props["Allocation"][mol]:
                data_alloc[mol][location+" - Count"] = props["Allocation"][mol][location][0]
                data_alloc[mol][location+" - mol/nm^2"] = form%props["Allocation"][mol][location][1]
                data_alloc[mol][location+" - mumol/m^2"] = form%props["Allocation"][mol][location][2]

        tables["alloc"] = pd.DataFrame.from_dict(data_alloc)

        # Time table
        tables["time"] = pd.DataFrame.from_dict(props["Generation_Time"], orient="index", columns={"Time (s)"}).round(decimals)

        # Return tables
        return tables


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

    def get_q(self):
        """Return the adjusted grid charges.

        Returns
        -------
        q : dictionary
            Dictionary of grid molecules and their charges
        """
        return self._q

    def get_drill(self):
        """Return the drilling axis.

        Returns
        -------
        drill : string
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
        props : dictionary
            Dictionary containing properties of the generated pore
        """
        return self._props
