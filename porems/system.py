################################################################################
# Basic Pore System Classes                                                    #
#                                                                              #
"""Here basic pore system constructions are defined."""
################################################################################


import porems as pms


class PoreSystem():
    """This class is a container class for further pore system classes."""
    def __init__(self):
        # Initialize
        self._sort_list = []
        self._site_in = []
        self._site_ex = []
        self._res = 0


    ################
    # Finalization #
    ################
    def finalize(self):
        """Finalize pore system."""
        # Fill silanol molecules on empty binding sites
        mols_in = self._pore.fill_sites(self._site_in, self._normal_in) if self._site_in else  []
        mols_ex = self._pore.fill_sites(self._site_ex, self._normal_ex) if self._site_ex else  []

        # Add to sorting list
        for mols in [mols_in, mols_ex]:
            for mol in mols:
                if not mol.get_short() in self._sort_list:
                    self._sort_list.append(mol.get_short())

        # Create reservoir
        self._pore.reservoir(self._res)

    def store(self, link="./"):
        """Store pore system and all necessary files for simulation at given
        link.

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


class PoreCylinder(PoreSystem):
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
    Following example generates a cylindrical pore with a diameter of 4nm,
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
        # Call super class
        super(PoreCylinder, self).__init__()

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


class PoreSlit(PoreSystem):
    """This class carves a slit-pore out of a :math:`\\beta`-cristobalite block.

    Parameters
    ----------
    size : list
        Size of the silicon-oxygen-grid
    height : float
        Pore height

    Examples
    --------
    Following example generates a slit-pore with a height of 4nm functionalized
    with TMS

    .. code-block:: python

        import porems as pms

        pore = pms.PoreSlit([6, 6, 6], 4)

        pore.attach(pms.gen.tms(), 0, [0, 1], 100, "in")

        pore.finalize()

        pore.store("output/")
    """
    def __init__(self, size, height):
        # Call super class
        super(PoreSlit, self).__init__()

        # Initialize
        self._size = size
        self._height = height
        self._sort_list = ["OM", "SI"]

        # Build pattern
        pattern = pms.BetaCristobalit()
        pattern.generate(self._size, "z")

        # Create block
        block = pattern.get_block()
        block.set_name("pore")

        # Dice up block
        dice = pms.Dice(block, 0.4, True)
        matrix = pms.Matrix(dice.find_parallel(None, ["Si", "O"], 0.155, 10e-2))

        # Carve out shape
        self._centroid = block.centroid()
        central = [0, 0, 1]
        self._cuboid = pms.Cuboid({"centroid": self._centroid, "central": central, "length": size[2], "width": size[0], "height": height})
        del_list = [atom_id for atom_id, atom in enumerate(block.get_atom_list()) if self._cuboid.is_in(atom.get_pos())]
        matrix.strip(del_list)

        # Prepare pore surface
        self._pore = pms.Pore(block, matrix)
        self._pore.set_name("pore")
        self._pore.prepare()

        # Determine sites
        self._pore.sites()
        site_list = self._pore.get_sites()
        self._site_in = [site_key for site_key, site_val in site_list.items() if site_val["type"]=="in"]

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
        return self._cuboid.normal(pos)

    def attach(self, mol, mount, axis, amount, scale=1, trials=1000, is_rotate=False):
        """Attach molecule on the surface.

        Parameters
        ----------
        mol : Molecule
            Molecule object to attach
        mount : integer
            Atom id of the molecule that is placed on the surface silicon atom
        axis : list
            List of two atom ids of the molecule that define the molecule axis
        scale : float, optional
            Circumference scaling around the molecule position
        trials : integer, optional
            Number of trials picking a random site
        is_rotate : bool, optional
            True to randomly rotate molecule around own axis
        """
        # Run attachment
        mols = self._pore.attach(mol, mount, axis, self._site_in, amount, self._normal_in, scale, trials, is_proxi=True, is_random=True, is_rotate=is_rotate)

        # Add to sorting list
        for mol in mols:
            if not mol.get_short() in self._sort_list:
                self._sort_list.append(mol.get_short())


class PoreCapsule(PoreSystem):
    """This class carves a capsule pore system out of a
    :math:`\\beta`-cristobalite block.

    Parameters
    ----------
    size : list
        Size of the silicon-oxygen-grid
    diam : float
        Cylinder diameter
    sep : float
        Seperation length between capsule
    res : float, optional
        Reservoir size on each side

    Examples
    --------
    Following example generates a capsule pore with a diameter of 4nm, a
    separation distance between the capsules of 2nm, reservoirs of 5nm on each
    side and a surface functionalized with TMS

    .. code-block:: python

        import porems as pms

        pore = pms.PoreCapsule([6, 6, 12], 4, 2)

        pore.attach(pms.gen.tms(), 0, [0, 1], 100, "in")
        pore.attach(pms.gen.tms(), 0, [0, 1], 20, "ex")

        pore.finalize()

        pore.store("output/")
    """
    def __init__(self, size, diam, sep, res=5):
        # Call super class
        super(PoreCapsule, self).__init__()

        # Initialize
        self._size = size
        self._diam = diam
        self._sep = sep
        self._res = res
        self._sort_list = ["OM", "SI"]

        # Build pattern
        self._pattern = pms.BetaCristobalit()
        self._pattern.generate(self._size, "z")
        self._pattern.exterior()
        self._len_cyl = (self._pattern.get_size()[2]-sep-diam)/2

        # Create block
        block = self._pattern.get_block()
        block.set_name("pore")

        # Dice up block
        dice = pms.Dice(block, 0.4, True)
        matrix = pms.Matrix(dice.find_parallel(None, ["Si", "O"], 0.155, 10e-2))
        oxygen_out = matrix.bound(1)

        # Carve out shape
        central = [0, 0, 1]

        self._centroid = {}
        self._centroid["block"] = block.centroid()
        self._centroid["cyl_l"] = self._centroid["block"][:2]+[0]
        self._centroid["cyl_r"] = self._centroid["block"][:2]+[self._pattern.get_size()[2]]
        self._centroid["sph_l"] = self._centroid["block"][:2]+[self._len_cyl]
        self._centroid["sph_r"] = self._centroid["block"][:2]+[self._pattern.get_size()[2]-self._len_cyl]

        self._shape = {}
        self._shape["cyl_l"] = pms.Cylinder({"centroid": self._centroid["cyl_l"], "central": central, "length": self._len_cyl*2, "diameter": diam})
        self._shape["cyl_r"] = pms.Cylinder({"centroid": self._centroid["cyl_r"], "central": central, "length": self._len_cyl*2, "diameter": diam})
        self._shape["sph_l"] = pms.Sphere({"centroid": self._centroid["sph_l"], "central": central, "diameter": diam})
        self._shape["sph_r"] = pms.Sphere({"centroid": self._centroid["sph_r"], "central": central, "diameter": diam})

        del_list = []
        for shape in self._shape.values():
            del_list.extend([atom_id for atom_id, atom in enumerate(block.get_atom_list()) if shape.is_in(atom.get_pos())])
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
        if pos[2] <= self._len_cyl:
            return self._shape["cyl_l"].normal(pos)
        elif pos[2] > self._len_cyl and pos[2] < self._centroid["block"][2]:
            return [x if i<2 else -x for i, x in enumerate(self._shape["sph_l"].normal(pos))]
        elif pos[2] > self._centroid["block"][2] and pos[2] < self._pattern.get_size()[2]-self._len_cyl:
            return [x if i<2 else -x for i, x in enumerate(self._shape["sph_r"].normal(pos))]
        elif pos[2] >= self._pattern.get_size()[2]-self._len_cyl:
            return self._shape["cyl_r"].normal(pos)

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
        return [0, 0, -1] if pos[2] < self._centroid["block"][2] else [0, 0, 1]

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
