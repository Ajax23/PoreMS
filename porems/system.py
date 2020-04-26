################################################################################
# Basic Pore System Classes                                                    #
#                                                                              #
"""Here basic pore system constructions are defined."""
################################################################################


import math
import pandas as pd
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
        mols_in = self._pore.fill_sites(self._site_in, self._normal_in, "in") if self._site_in else  []
        mols_ex = self._pore.fill_sites(self._site_ex, self._normal_ex, "ex") if self._site_ex else  []

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
        self._diam = diam
        self._res = res
        self._sort_list = ["OM", "SI"]

        # Build pattern
        pattern = pms.BetaCristobalit()
        pattern.generate(size, "z")
        pattern.exterior()

        # Create block
        self._block = pattern.get_block()
        self._block.set_name("pore")
        self._box = self._block.get_box()

        # Dice up block
        dice = pms.Dice(self._block, 0.4, True)
        matrix = pms.Matrix(dice.find_parallel(None, ["Si", "O"], 0.155, 10e-2))
        oxygen_out = matrix.bound(1)

        # Carve out shape
        self._centroid = self._block.centroid()
        central = [0, 0, 1]
        self._cylinder = pms.Cylinder({"centroid": self._centroid, "central": central, "length": size[2], "diameter": diam-0.5})  # Preperation precaution
        del_list = [atom_id for atom_id, atom in enumerate(self._block.get_atom_list()) if self._cylinder.is_in(atom.get_pos())]
        matrix.strip(del_list)

        # Prepare pore surface
        self._pore = pms.Pore(self._block, matrix)
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


    ##############
    # Attachment #
    ##############
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

    def attach(self, mol, mount, axis, amount, site_type, scale=1, trials=1000, inp="num", is_rotate=False):
        """Attach molecule on the surface.

        Parameters
        ----------
        mol : Molecule
            Molecule object to attach
        mount : integer
            Atom id of the molecule that is placed on the surface silicon atom
        axis : list
            List of two atom ids of the molecule that define the molecule axis
        amount : int
            Number of molecules to attach
        site_type : string
            Use **in** for the interior surface and **ex** for the exterior
        scale : float, optional
            Circumference scaling around the molecule position
        trials : integer, optional
            Number of trials picking a random site
        inp : string, optional
            Input type - **num** Number of molecules, **molar** :math:`\\frac{\\mu\\text{mol}}{\\text{m}^2}`
        is_rotate : bool, optional
            True to randomly rotate molecule around own axis
        """
        # Process input
        if site_type not in ["in", "ex"]:
            print("Pore: Wrong site_type input...")
            return

        if inp not in ["num", "molar"]:
            print("Pore: Wrong inp type...")
            return

        # Amount
        amount = int(pms.utils.mumol_m2_to_mols(amount, self.surface()[site_type])) if inp=="molar" else amount

        # Sites and normal vector
        sites = self._site_in if site_type=="in" else self._site_ex
        normal = self._normal_in if site_type=="in" else self._normal_ex

        # Run attachment
        mols = self._pore.attach(mol, mount, axis, sites, amount, normal, scale, trials, site_type=site_type, is_proxi=True, is_random=True, is_rotate=is_rotate)

        # Add to sorting list
        for mol in mols:
            if not mol.get_short() in self._sort_list:
                self._sort_list.append(mol.get_short())

    def attach_special(self, mol, mount, axis, amount, scale=1, symmetry="point", is_rotate=False):
        """Special attachment of molecules on the surface.

        Parameters
        ----------
        mol : Molecule
            Molecule object to attach
        mount : integer
            Atom id of the molecule that is placed on the surface silicon atom
        axis : list
            List of two atom ids of the molecule that define the molecule axis
        amount : int
            Number of molecules to attach
        scale : float, optional
            Circumference scaling around the molecule position
        symmetry : string, optional
            Symmetry option - point, mirror
        is_rotate : bool, optional
            True to randomly rotate molecule around own axis
        """
        # Process input
        if symmetry not in ["point", "mirror"]:
            print("Symmetry type not supported...")
            return

        # Calculate geometrical positions
        dist = self._box[2]/amount if amount>0 else 0
        start = dist/2

        pos_list = []
        for i in range(amount):
            if symmetry == "point":
                coeff = -1 if i % 2 == 0 else 1
            elif symmetry == "mirror":
                coeff = 1

            x = self._centroid[0]+coeff*self.diameter()/2
            y = self._centroid[1]
            z = start+dist*i

            pos_list.append([x, y, z])

        # Run attachment
        mols = self._pore.attach(mol, mount, axis, self._site_in, len(pos_list), self._normal_in, scale, pos_list=pos_list, is_proxi=True, is_random=False, is_rotate=is_rotate)

        # Add to sorting list
        for mol in mols:
            if not mol.get_short() in self._sort_list:
                self._sort_list.append(mol.get_short())


    ############
    # Analysis #
    ############
    def diameter(self):
        """Calculate true cylinder diameter after drilling and preperation. This
        is done by determining the mean value :math:`\\bar r` of the silicon
        distances :math:`r_i` of silicon :math:`i` towards the pore center

        .. math::

            \\bar r=\\frac1n\\sum_{i=1}^nr_i

        with the number of silicon atoms :math:`n`. The diameter is then

        .. math::

            d=2\\bar r=\\frac2n\\sum_{i=1}^nr_i.

        Returns
        -------
        diameter : float
            Pore diameter after preperation
        """
        # Calculate distance towards central axis of binding site silicon atoms
        r = [pms.geom.length(pms.geom.vector([self._centroid[0], self._centroid[1], self._block.pos(site)[2]], self._block.pos(site))) for site in self._site_in]

        # Calculate mean
        r_bar = sum(r)/len(r) if len(r)>0 else 0

        # Calculate diameter
        return 2*r_bar

    def roughness(self):
        """Calculate surface roughness. In the case of a cylindric pore one can
        visualize pulling the pore apart, thus flattening the interior surface.
        The roughness is then determined by calculating the standard deviation
        of the binding site silicon atoms peaks and valleys.

        It is therefore sufficient to calculate the distnaces towards a specific
        axis, whis in this case will be the central axis. The mean value
        :math:`\\bar r` of the silicon distances :math:`r_i` of silicon
        :math:`i` towards the pore centre, is calculated by

        .. math::

            \\bar r=\\frac1n\\sum_{i=1}^nr_i

        with the number of silicon atoms :math:`n`. This mean value is used in
        the square root roughness calculation

        .. math::

            R_q = \\sqrt{\\frac1n\\sum_{i=1}^n\\|r_i-\\bar r\\|^2}.

        Returns
        -------
        roughness : float
            Surface roughness
        """
        # Calculate distance towards central axis of binding site silicon atoms
        r = [pms.geom.length(pms.geom.vector([self._centroid[0], self._centroid[1], self._block.pos(site)[2]], self._block.pos(site))) for site in self._site_in]

        # Calculate mean
        r_bar = sum(r)/len(r) if len(r)>0 else 0

        # Calculate square root roughness
        return math.sqrt(sum([(r_i-r_bar)**2 for r_i in r])/len(r)) if len(r)>0 else 0

    def volume(self):
        """Calculate pore volume. This is done by defining a new shape object
        with system sizes after pore preperation and using the volume function
        :func:`porems.shape.Cylinder.volume`.

        Returns
        -------
        volume : float
            Pore volume
        """
        return pms.Cylinder({"centroid": self._centroid, "central": [0, 0, 1], "length": self._box[2], "diameter": self.diameter()}).volume()

    def surface(self):
        """Calculate pore surface and exterior surface. This is done by defining
        a new shape object with system sizes after pore preperation and using
        the surface function :func:`porems.shape.Cylinder.surface`.

        Returns
        -------
        surface : dictionary
            Pore surface of interior and exterior
        """
        diam = self.diameter()

        surf_in = pms.Cylinder({"centroid": self._centroid, "central": [0, 0, 1], "length": self._box[2], "diameter": diam}).surface()
        surf_ex = 2*(self._box[0]*self._box[1]-math.pi*(diam/2)**2)

        return {"in": surf_in, "ex": surf_ex}

    def allocation(self):
        """Calculate molecule allocation on the surface. Using interior and
        exterior surfaces, the allocation rates can be determined by counting
        the number of used molecules on the surfaces.

        Using the conversion function :func:`porems.utils.mols_to_mumol_m2`, the
        number of molecules is converted to concentration in
        :math:`\\frac{\\mu\\text{mol}}{\\text{m}^2}`.

        Returns
        -------
        alloc : dictionary
            Dictionary containing the surface allocation of all molecules in
            number of molecules, :math:`\\frac{\\text{mols}}{\\text{nm}^2}` and
            :math:`\\frac{\\mu\\text{mol}}{\\text{m}^2}`
        """
        # Get surfaces
        surf = self.surface()
        site_dict = self._pore.get_site_dict()

        # Calculate allocation for all molecules
        alloc = {}
        for mol in sorted(self._sort_list):
            for site_type in ["in", "ex"]:
                if mol in site_dict[site_type]:
                    if not mol in alloc:
                        alloc[mol] = {"in": [0, 0, 0], "ex": [0, 0, 0]}
                    # Number of molecules
                    alloc[mol][site_type][0] = len(site_dict[site_type][mol])

                    # Molecules per nano meter
                    alloc[mol][site_type][1] = len(site_dict[site_type][mol])/surf[site_type] if surf[site_type]>0 else 0

                    # Micromolar per meter
                    alloc[mol][site_type][2] = pms.utils.mols_to_mumol_m2(len(site_dict[site_type][mol]), surf[site_type]) if surf[site_type]>0 else 0

        # OH allocation
        alloc["OH"] = {"in": [0, 0, 0], "ex": [0, 0, 0]}
        for site_type in ["in", "ex"]:
            num_oh = 0
            num_oh += len(site_dict[site_type]["SL"]) if "SL" in site_dict[site_type] else 0
            num_oh += len(site_dict[site_type]["SLG"])*2 if "SLG" in site_dict[site_type] else 0

            alloc["OH"][site_type][0] = num_oh
            alloc["OH"][site_type][1] = num_oh/surf[site_type] if surf[site_type]>0 else 0
            alloc["OH"][site_type][2] = pms.utils.mols_to_mumol_m2(num_oh, surf[site_type]) if surf[site_type]>0 else 0

        # Hydroxilation - Total number of binding sites
        alloc["Hydro"] = {"in": [0, 0, 0], "ex": [0, 0, 0]}
        for site_type in ["in", "ex"]:
            num_tot = len(sum([x["o"] for x in self._pore.get_sites().values() if x["type"]==site_type], []))
            alloc["Hydro"][site_type][0] = num_tot
            alloc["Hydro"][site_type][1] = num_tot/surf[site_type] if surf[site_type]>0 else 0
            alloc["Hydro"][site_type][2] = pms.utils.mols_to_mumol_m2(num_tot, surf[site_type]) if surf[site_type]>0 else 0

        return alloc

    def table(self, decimals=3):
        """Create properties as pandas table for easy viewing.

        Parameters
        ----------
        decimals : integer, optional
            Number of decimals to be rounded to

        Returns
        -------
        tables : dictionary
            Dictionary of pandas table of all properties
        """
        # Initialize
        tables = {}
        form = "%."+str(decimals)+"f"

        # Properties table
        data_props = {}
        data_props["Dimension"] = "["+form%self._block.get_box()[0]+", "+form%self._block.get_box()[1]+", "+form%self._block.get_box()[2]+"]"
        data_props["Diameter"] = form%self.diameter()
        data_props["Reservoir"] = form%self._res
        data_props["Roughness"] = form%self.roughness()
        data_props["Surface"] = form%self.surface()["in"]
        data_props["Volume"] = form%self.volume()

        tables["props"] = pd.DataFrame.from_dict(data_props, orient="index", columns={"d = "+"%.2f"%self._diam+" nm"})

        # Allocation table
        data_alloc = {}
        allocation = self.allocation()
        for mol in allocation:
            data_alloc[mol] = {}
            for site_type in allocation[mol]:
                data_alloc[mol][site_type+" - Count"] = allocation[mol][site_type][0]
                data_alloc[mol][site_type+" - mols/nm^2"] = form%allocation[mol][site_type][1]
                data_alloc[mol][site_type+" - mumol/m^2"] = form%allocation[mol][site_type][2]

        tables["alloc"] = pd.DataFrame.from_dict(data_alloc)

        return tables


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


    ##############
    # Attachment #
    ##############
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
        amount : int
            Number of molecules to attach
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
        Separation length between capsule
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


    ##############
    # Attachment #
    ##############
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
        amount : int
            Number of molecules to attach
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
        mols = self._pore.attach(mol, mount, axis, sites, amount, normal, scale, trials, site_type=site_type, is_proxi=True, is_random=True, is_rotate=is_rotate)

        # Add to sorting list
        for mol in mols:
            if not mol.get_short() in self._sort_list:
                self._sort_list.append(mol.get_short())
