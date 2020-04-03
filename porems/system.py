################################################################################
# Basic Pore System Classes                                                    #
#                                                                              #
"""Here basic pore system constructions are defined."""
################################################################################


import porems as pms


class PoreCylinder():
    """This class carves a cylindric pore system out of a
    :math:`\\beta`-cristobalite block.

    Parameters
    ----------
    size : list
        Size of the silicon-oxygen-grid
    diam : float
        Cylinder diameter
    res : float, optional
        Reservoir size on each side

    Examples
    --------
    Following example generates a cylindrical pore with a diameter of 4, a
    reservoirs of 5nm on each side and a surface functionalized with TMS

    .. code-block:: python

        import porems as pms

        pore = pms.PoreCylinder([6, 6, 6], 4, 5)

        pore.attach(pms.gen.tms(), 0, [0, 1], 100, "in")
        pore.attach(pms.gen.tms(), 0, [0, 1], 20, "ex")

        pore.finalize()

        pore.store("output/")
    """
    def __init__(self, size, diam, res=5):
        # Initialize
        self._size = size
        self._diam = diam
        self._res = res
        self._sort_list = ["OM", "SI"]

        # Build pattern
        pattern = pms.BetaCristobalit()
        pattern.generate(self._size, "z")
        pattern.exterior()

        # Create block
        block = pattern.get_block()
        block.set_name("pore")

        # Dice up block
        dice = pms.Dice(block, 0.4, True)
        matrix = pms.Matrix(dice.find_parallel(None, ["Si", "O"], 0.155, 10e-2))
        oxygen_out = matrix.bound(1)

        # Carve out shape
        self._centroid = block.centroid()
        central = [0, 0, 1]
        self._cylinder = pms.Cylinder({"centroid": self._centroid, "central": central, "length": size[2], "diameter": diam})
        del_list = [atom_id for atom_id, atom in enumerate(block.get_atom_list()) if self._cylinder.is_in(atom.get_pos())]
        matrix.strip(del_list)

        # Prepare pore surface
        self._pore = pms.Pore(block, matrix)
        self._pore.set_name("pore")
        self._pore.prepare()

        # Determine sites
        self._pore.sites(oxygen_out)
        site_list = self._pore.get_sites()
        self._site_in = [site_key for site_key, site_val in site_list.items() if site_val["type"]=="in"]
        self._site_ex = [site_key for site_key, site_val in site_list.items() if site_val["type"]=="ex"]

        # Objectify grid
        grid_atoms = [atom for atom in matrix.bound(0, "gt") if not atom in matrix.bound(1)+list(site_list.keys())]
        self._pore.objectify(grid_atoms)


    ###############
    # Attachement #
    ###############
    def _normal_in(self, pos):
        """Normal function for the interior surface

        Parameters
        ----------
        pos : list
            Position on the surface

        Returns
        -------
        normal : list
            Vector perpendicular to surface of the given position
        """
        return self._cylinder.normal(pos)

    def _normal_ex(self, pos):
        """Normal function for the exterior surface

        Parameters
        ----------
        pos : list
            Position on the surface

        Returns
        -------
        normal : list
            Vector perpendicular to surface of the given position
        """
        return [0, 0, -1] if pos[2] < self._centroid[2] else [0, 0, 1]

    def attach(self, mol, mount, axis, amount, site_type, scale=1, trials=1000, is_rotate=False):
        """Attach molecule on the surface.

        Parameters
        ----------
        mol : Molecule
            Molecule object to attach
        mount : integer
            Atom id of the molecule that is placed on the surface silicon atom
        axis : list
            List of two atom ids of the molecule that define the molecule axis
        site_type : string
            Use **in** for the interior surface and **ex** for the exterior
        scale : float, optional
            Circumference scaling around the molecule position
        trials : integer, optional
            Number of trials picking a random site
        is_rotate : bool, optional
            True to randomly rotate molecule around own axis
        """
        # Process input
        if site_type not in ["in", "ex"]:
            print("Pore: Wrong site_type input...")
            return

        sites = self._site_in if site_type=="in" else self._site_ex
        normal = self._normal_in if site_type=="in" else self._normal_ex

        # Run attachment
        mols = self._pore.attach(mol, mount, axis, sites, amount, normal, scale, trials, is_proxi=True, is_random=True, is_rotate=is_rotate)

        # Add to sorting list
        for mol in mols:
            if not mol.get_short() in self._sort_list:
                self._sort_list.append(mol.get_short())


    ################
    # Finalization #
    ################
    def finalize(self):
        """Finalize pore system."""
        # Fill silanol molecules on empty binding sites
        mols_in = self._pore.fill_sites(self._site_in, self._normal_in)
        mols_ex = self._pore.fill_sites(self._site_ex, self._normal_ex)

        # Add to sorting list
        for mols in [mols_in, mols_ex]:
            for mol in mols:
                if not mol.get_short() in self._sort_list:
                    self._sort_list.append(mol.get_short())

        # Create reservoir
        self._pore.reservoir(self._res)

    def store(self, link="./"):
        """Store pore system and all necessary files for simulation at given link.

        Parameters
        ----------
        link : string, optional
            Folder link for output
        """
        # Process input
        self._link = link if link[-1] == "/" else link+"/"

        # Create store object
        store = pms.Store(self._pore, link, self._sort_list)

        # Save files
        store.gro(use_atom_names=True)
        store.obj()
        store.top()
        store.grid()
