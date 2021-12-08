################################################################################
# Basic Pore System Classes                                                    #
#                                                                              #
"""Here basic pore system constructions are defined."""
################################################################################


import os
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
        self._pore_shape = ""


    ##############
    # Attachment #
    ##############
    def _siloxane(self, hydro, site_type, slx_dist=[0.507-1e-2, 0.507+1e-2]):
        """Attach siloxane bridges using function
        :func:`porems.pore.Pore.siloxane`.

        Parameters
        ----------
        hydro: float
            Hydroxilation degree in
            :math:`\\frac{\\mu\\text{mol}}{\\text{m}^2}`.
        slx_dist : list
            Silicon atom distance to search for parters in proximity
        site_type : string, optional
            Site type - interior **in**, exterior **ex**
        """
        # Initialize
        site_list = self._pore.get_sites()
        sites = self._site_in if site_type=="in" else self._site_ex
        normal = self._normal_in if site_type=="in" else self._normal_ex

        # Amount - Connect two oxygen to one siloxane
        oh = len(sum([site_list[site]["o"] for site in sites], []))
        oh_goal = pms.utils.mumol_m2_to_mols(hydro, self.surface()[site_type])
        amount = round((oh-oh_goal)/2)

        # Fill siloxane
        if amount > 0:
            # Run attachment
            mols = self._pore.siloxane(sites, amount, normal, slx_dist=slx_dist, site_type=site_type)

            # Add to sorting list
            for mol in mols:
                if not mol.get_short() in self._sort_list:
                    self._sort_list.append(mol.get_short())

    def attach(self, mol, mount, axis, amount, site_type="in", inp="num", pos_list=[], scale=1, trials=1000, is_proxi=True, is_rotate=False):
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
        site_type : string, optional
            Use **in** for the interior surface and **ex** for the exterior
        inp : string, optional
            Input type: **num** - Number of molecules,
            **molar** - :math:`\\frac{\\mu\\text{mol}}{\\text{m}^2}`,
            **percent** - :math:`\\%` of OH groups
        pos_list : list, optional
            List of positions (cartesian) to find nearest available binding site for
        scale : float, optional
            Circumference scaling around the molecule position
        trials : integer, optional
            Number of trials picking a random site
        is_proxi : bool, optional
            True to fill binding sites in proximity of filled binding site
        is_rotate : bool, optional
            True to randomly rotate molecule around own axis
        """
        # Process input
        if site_type not in ["in", "ex"]:
            print("Pore: Wrong site_type input...")
            return

        if inp not in ["num", "molar", "percent"]:
            print("Pore: Wrong inp type...")
            return

        # Sites and normal vector
        sites = self._site_in if site_type=="in" else self._site_ex
        normal = self._normal_in if site_type=="in" else self._normal_ex

        # Amount
        if inp=="molar":
            amount = int(pms.utils.mumol_m2_to_mols(amount, self.surface()[site_type]))
        elif inp=="percent":
            num_oh = len(sites)
            num_oh += sum([1 for x in self._pore.get_sites().values() if len(x["o"])==2 and x["type"]==site_type])
            amount = int(amount/100*num_oh)

        # Check number of given positions
        if pos_list and not len(pos_list)==amount:
            print("Pore: Number of given positions does not match number of groups to attach...")
            return

        # Run attachment
        mols = self._pore.attach(mol, mount, axis, sites, amount, normal, scale, trials, pos_list=pos_list, site_type=site_type, is_proxi=is_proxi, is_random=True, is_rotate=is_rotate)

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
        mols_in = self._pore.fill_sites(self._site_in, self._normal_in, "in") if self._site_in else []
        mols_ex = self._pore.fill_sites(self._site_ex, self._normal_ex, "ex") if self._site_ex else []

        # Add to sorting list
        for mols in [mols_in, mols_ex]:
            for mol in mols:
                if not mol.get_short() in self._sort_list:
                    self._sort_list.append(mol.get_short())

        # Create reservoir
        self._pore.reservoir(self._res)

    def store(self, link="./", sort_list=[]):
        """Store pore system and all necessary files for simulation at given
        link.

        Parameters
        ----------
        link : string, optional
            Folder link for output
        """
        # Process input
        link = link if link[-1] == "/" else link+"/"

        # Set sort list
        sort_list = sort_list if sort_list else self._sort_list

        # Create store object
        store = pms.Store(self._pore, link, sort_list)

        # Save files
        store.gro(use_atom_names=True)
        store.obj()
        store.top()
        store.grid()
        pms.utils.save(self, link+self._pore.get_name()+"_system.obj")


    ############
    # Analysis #
    ############
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
        # num_in_ex = self._pore.get_num_in_ex()

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
            num_oh = len(sum([x["o"] for x in self._pore.get_sites().values() if x["type"]==site_type], []))
            for mol in site_dict[site_type].keys():
                num_oh -= len(site_dict[site_type][mol]) if mol not in ["SL", "SLG", "SLX"] else 0

            # num_oh = num_oh-num_in_ex if site_type=="ex" else num_oh+num_in_ex

            alloc["OH"][site_type][0] = num_oh
            alloc["OH"][site_type][1] = num_oh/surf[site_type] if surf[site_type]>0 else 0
            alloc["OH"][site_type][2] = pms.utils.mols_to_mumol_m2(num_oh, surf[site_type]) if surf[site_type]>0 else 0

        # Hydroxylation - Total number of binding sites
        alloc["Hydro"] = {"in": [0, 0, 0], "ex": [0, 0, 0]}
        for site_type in ["in", "ex"]:
            num_tot = len(sum([x["o"] for x in self._pore.get_sites().values() if x["type"]==site_type], []))

            # num_tot = num_tot-num_in_ex if site_type=="ex" else num_tot+num_in_ex

            alloc["Hydro"][site_type][0] = num_tot
            alloc["Hydro"][site_type][1] = num_tot/surf[site_type] if surf[site_type]>0 else 0
            alloc["Hydro"][site_type][2] = pms.utils.mols_to_mumol_m2(num_tot, surf[site_type]) if surf[site_type]>0 else 0

        return alloc

    def reservoir(self):
        """Return the reservoir length.

        Returns
        -------
        res : float
            Reservoir length
        """
        return self._res

    def box(self):
        """Return the box size of the pore block.

        Returns
        -------
        box : list
            Box size in all dimensions
        """
        return self._pore.get_box()

    def centroid(self):
        """Return pore centroid.

        Returns
        -------
        centroid : list
            Centroid of the pore
        """
        return self._centroid

    def shape(self):
        """Return the pore shape for analysis using PoreAna.

        Returns
        -------
        pore_shape : string
            Pore shape
        """
        return self._pore_shape


    #########
    # Table #
    #########
    def _table_base(self, props, decimals=3):
        """Base functino for converting properties to pandas table for easy viewing.

        Parameters
        ----------
        props : dictionary
            Additional pore specific properties {propname: [in-prop, ex-prop], ...}
        decimals : integer, optional
            Number of decimals to be rounded to

        Returns
        -------
        tables : DataFrame
            Pandas table of all properties
        """
        # Initialize
        form = "%."+str(decimals)+"f"

        # Get allocation data
        allocation = self.allocation()

        # Save data
        data = {"Interior": {}, "Exterior": {}}

        data["Interior"]["Silica block xyz-dimensions (nm)"] = " "
        data["Exterior"]["Silica block xyz-dimensions (nm)"] = "["+form%self.box()[0]+", "+form%self.box()[1]+", "+form%(self.box()[2]-2*self.reservoir())+"]"
        data["Interior"]["Simulation box xyz-dimensions (nm)"] = " "
        data["Exterior"]["Simulation box xyz-dimensions (nm)"] = "["+form%self.box()[0]+", "+form%self.box()[1]+", "+form%self.box()[2]+"]"
        data["Interior"]["Pore drilling direction"] = "z"
        data["Exterior"]["Pore drilling direction"] = " "
        data["Interior"]["Surface roughness (nm)"] = form%self.roughness()["in"] if "in" in self.roughness() else form%0
        data["Exterior"]["Surface roughness (nm)"] = form%self.roughness()["ex"] if "ex" in self.roughness() else form%0

        for prop_name, values in props.items():
            data["Interior"][prop_name] = values[0]
            data["Exterior"][prop_name] = values[1]

        data["Interior"]["Solvent reservoir z-dimension (nm)"] = " "
        data["Exterior"]["Solvent reservoir z-dimension (nm)"] = form%self.reservoir()
        data["Interior"]["Pore volume (nm^3)"] = form%self.volume()
        data["Exterior"]["Pore volume (nm^3)"] = " "
        data["Interior"]["Solvent reservoir volume (nm^3)"] = " "
        data["Exterior"]["Solvent reservoir volume (nm^3)"] = "2 * "+form%(self.box()[0]*self.box()[1]*self.reservoir())
        data["Interior"]["Surface area (nm^2)"] = form%self.surface()["in"]
        data["Exterior"]["Surface area (nm^2)"] = "2 * "+form%(self.surface()["ex"]/2)

        data["Interior"]["Surface chemistry - Before Functionalization"] = " "
        data["Exterior"]["Surface chemistry - Before Functionalization"] = " "
        data["Interior"]["    Number of single silanol groups"] = "%i"%sum([1 for x in self._pore.get_sites().values() if len(x["o"])==1 and x["type"]=="in"])
        data["Exterior"]["    Number of single silanol groups"] = "%i"%sum([1 for x in self._pore.get_sites().values() if len(x["o"])==1 and x["type"]=="ex"])
        data["Interior"]["    Number of geminal silanol groups"] = "%i"%sum([1 for x in self._pore.get_sites().values() if len(x["o"])==2 and x["type"]=="in"])
        data["Exterior"]["    Number of geminal silanol groups"] = "%i"%sum([1 for x in self._pore.get_sites().values() if len(x["o"])==2 and x["type"]=="ex"])
        data["Interior"]["    Number of siloxane bridges"] = "%i"%allocation["SLX"]["in"][0] if "SLX" in allocation else "0"
        data["Exterior"]["    Number of siloxane bridges"] = "%i"%allocation["SLX"]["ex"][0] if "SLX" in allocation else "0"
        data["Interior"]["    Total number of OH groups"] = "%i"%allocation["Hydro"]["in"][0]
        data["Exterior"]["    Total number of OH groups"] = "%i"%allocation["Hydro"]["ex"][0]
        data["Interior"]["    Overall hydroxylation (mumol/m^2)"] = form%allocation["Hydro"]["in"][2]
        data["Exterior"]["    Overall hydroxylation (mumol/m^2)"] = form%allocation["Hydro"]["ex"][2]

        data["Interior"]["Surface chemistry - After Functionalization"] = " "
        data["Exterior"]["Surface chemistry - After Functionalization"] = " "
        for mol in allocation.keys():
            if mol not in ["SL", "SLG", "SLX", "Hydro", "OH"]:
                data["Interior"]["    Number of "+mol+" groups"] = "%i"%allocation[mol]["in"][0]
                data["Exterior"]["    Number of "+mol+" groups"] = "%i"%allocation[mol]["ex"][0]
                data["Interior"]["    "+mol+" density (mumol/m^2)"] = form%allocation[mol]["in"][2]
                data["Exterior"]["    "+mol+" density (mumol/m^2)"] = form%allocation[mol]["ex"][2]
        data["Interior"]["    Bonded-phase density (mumol/m^2)"] = form%(allocation["Hydro"]["in"][2]-allocation["OH"]["in"][2])
        data["Exterior"]["    Bonded-phase density (mumol/m^2)"] = form%(allocation["Hydro"]["ex"][2]-allocation["OH"]["ex"][2])
        data["Interior"]["    Number of residual OH groups"] = "%i"%allocation["OH"]["in"][0]
        data["Exterior"]["    Number of residual OH groups"] = "%i"%allocation["OH"]["ex"][0]
        data["Interior"]["    Residual hydroxylation (mumol/m^2)"] = form%allocation["OH"]["in"][2]
        data["Exterior"]["    Residual hydroxylation (mumol/m^2)"] = form%allocation["OH"]["ex"][2]

        return pd.DataFrame.from_dict(data)


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
    hydro: list, optional
        Hydroxilation degree for interior and exterior of the pore in
        :math:`\\frac{\\mu\\text{mol}}{\\text{m}^2}`

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
    def __init__(self, size, diam, res=5, hydro=[0, 0]):
        # Call super class
        super(PoreCylinder, self).__init__()

        # Initialize
        self._diam = diam
        self._res = res
        self._sort_list = ["OM", "SI"]

        # Build pattern
        pattern = pms.BetaCristobalit()
        pattern.generate(size, "z")

        # Create block
        self._block = pattern.get_block()
        self._block.set_name("pore")
        self._box = self._block.get_box()

        # Dice up block
        dice = pms.Dice(self._block, 0.4, True)
        matrix = pms.Matrix(dice.find_parallel(None, ["Si", "O"], [0.155-1e-2, 0.155+1e-2]))

        # Create pore object
        self._pore = pms.Pore(self._block, matrix)
        self._pore.set_name("pore")

        # Create exterior
        if res:
            self._pore.exterior(pattern.get_gap())

        # Carve out shape
        self._centroid = self._block.centroid()
        central = [0, 0, 1]
        self._cylinder = pms.Cylinder({"centroid": self._centroid, "central": central, "length": size[2], "diameter": diam-0.5})  # Preperation precaution
        del_list = [atom_id for atom_id, atom in enumerate(self._block.get_atom_list()) if self._cylinder.is_in(atom.get_pos())]
        matrix.strip(del_list)
        self._pore_shape = "CYLINDER"

        # Prepare pore surface
        self._pore.prepare()

        # Determine sites
        self._pore.sites()
        site_list = self._pore.get_sites()
        self._site_in = [site_key for site_key, site_val in site_list.items() if site_val["type"]=="in"]
        self._site_ex = [site_key for site_key, site_val in site_list.items() if site_val["type"]=="ex"]

        # Save original surface si positions
        self._si_pos_in = [self._block.pos(site_key) for site_key, site_val in site_list.items() if site_val["type"]=="in"]
        self._si_pos_ex = [self._block.pos(site_key) for site_key, site_val in site_list.items() if site_val["type"]=="ex"]

        # Siloxane bridges
        if hydro[0]:
            self._siloxane(hydro[0], "in")
            self._site_in = [site_key for site_key, site_val in site_list.items() if site_val["type"]=="in"]
        if hydro[1]:
            self._siloxane(hydro[1], "ex")
            self._site_ex = [site_key for site_key, site_val in site_list.items() if site_val["type"]=="ex"]

        # Objectify grid
        non_grid = matrix.bound(1)+list(site_list.keys())
        bonded = matrix.bound(0, "gt")
        grid_atoms = [atom for atom in bonded if not atom in non_grid]
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

    def attach_special(self, mol, mount, axis, amount, scale=1, symmetry="point", is_proxi=True, is_rotate=False):
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
        is_proxi : bool, optional
            True to fill binding sites in proximity of filled binding site
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
        mols = self._pore.attach(mol, mount, axis, self._site_in, len(pos_list), self._normal_in, scale, pos_list=pos_list, is_proxi=is_proxi, is_random=False, is_rotate=is_rotate)

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
        r = [pms.geom.length(pms.geom.vector([self._centroid[0], self._centroid[1], pos[2]], pos)) for pos in self._si_pos_in]

        # Calculate mean
        r_bar = sum(r)/len(r) if len(r)>0 else 0

        # Calculate diameter
        return 2*r_bar

    def roughness(self):
        """Calculate surface roughness. In the case of a cylindric pore one can
        visualize pulling the pore apart, thus flattening the interior surface.
        The roughness is then determined by calculating the standard deviation
        of the binding site silicon atoms peaks and valleys.

        It is therefore enough to calculate the distances towards a specific
        axis, which in this case will be the central axis. The mean value
        :math:`\\bar r` of the silicon distances :math:`r_i` of silicon
        :math:`i` towards the pore center, is calculated by

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
        # Interior
        ## Calculate distance towards central axis of binding site silicon atoms
        r_in = [pms.geom.length(pms.geom.vector([self._centroid[0], self._centroid[1], pos[2]], pos)) for pos in self._si_pos_in]

        # Exterior
        if self._res:
            ## Create molecules with exterior positions
            temp_mol = pms.Molecule()
            for pos in self._si_pos_ex:
                temp_mol.add("Si", pos)
            temp_mol.zero()
            size = temp_mol.get_box()[2]

        ## Calculate distance to boundary
        r_ex = [pos[2] if pos[2] < size/2 else abs(pos[2]-size) for pos in self._si_pos_ex]

        # Calculate mean
        r_bar_in = sum(r_in)/len(r_in) if len(r_in)>0 else 0
        r_bar_ex = sum(r_ex)/len(r_ex) if len(r_ex)>0 else 0

        # Calculate roughness
        r_q_in =  math.sqrt(sum([(r_i-r_bar_in)**2 for r_i in r_in])/len(r_in)) if len(r_in)>0 else 0
        r_q_ex =  math.sqrt(sum([(r_i-r_bar_ex)**2 for r_i in r_ex])/len(r_ex)) if len(r_ex)>0 else 0

        # Calculate square root roughness
        return {"in": r_q_in, "ex": r_q_ex}

    def volume(self):
        """Calculate pore volume. This is done by defining a new shape object
        with system sizes after pore preparation and using the volume function
        :func:`porems.shape.Cylinder.volume`.

        Returns
        -------
        volume : float
            Pore volume
        """
        return pms.Cylinder({"centroid": self._centroid, "central": [0, 0, 1], "length": self._box[2], "diameter": self.diameter()}).volume()

    def surface(self):
        """Calculate pore surface and exterior surface. This is done by defining
        a new shape object with system sizes after pore preparation and using
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

    def table(self, decimals=3):
        """Create properties as pandas table for easy viewing.

        Parameters
        ----------
        decimals : integer, optional
            Number of decimals to be rounded to

        Returns
        -------
        tables : DataFrame
            Pandas table of all properties
        """
        # Initialize
        form = "%."+str(decimals)+"f"

        # Define pore specific properties
        props = {}
        props["Pore diameter (nm)"] = [form%self.diameter(), " "]

        return self._table_base(props, decimals)


class PoreSlit(PoreSystem):
    """This class carves a slit-pore out of a :math:`\\beta`-cristobalite block.

    Parameters
    ----------
    size : list
        Size of the silicon-oxygen-grid
    height : float
        Pore height
    res : float, optional
        Reservoir size on each side
    hydro: list, optional
        Hydroxilation degree for interior and exterior of the pore in
        :math:`\\frac{\\mu\\text{mol}}{\\text{m}^2}`

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
    def __init__(self, size, height, res=0, hydro=[0, 0]):
        # Call super class
        super(PoreSlit, self).__init__()

        # Initialize
        self._size = size
        self._height = height
        self._res = res
        self._sort_list = ["OM", "SI"]

        # Build pattern
        pattern = pms.BetaCristobalit()
        pattern.generate(size, "z")

        # Create block
        self._block = pattern.get_block()
        self._block.set_name("pore")
        self._box = self._block.get_box()

        # Dice up block
        dice = pms.Dice(self._block, 0.4, True)
        matrix = pms.Matrix(dice.find_parallel(None, ["Si", "O"], [0.155-1e-2, 0.155+1e-2]))

        # Create pore object
        self._pore = pms.Pore(self._block, matrix)
        self._pore.set_name("pore")

        # Create exterior
        if res:
            self._pore.exterior(pattern.get_gap())

        # Carve out shape
        self._centroid = self._block.centroid()
        central = [0, 0, 1]
        self._cuboid = pms.Cuboid({"centroid": self._centroid, "central": central, "length": size[2], "width": size[0], "height": height-0.5})  # Preperation precaution
        del_list = [atom_id for atom_id, atom in enumerate(self._block.get_atom_list()) if self._cuboid.is_in(atom.get_pos())]
        matrix.strip(del_list)
        self._pore_shape = "SLIT"

        # Prepare pore surface
        self._pore.prepare()

        # Determine sites
        self._pore.sites()
        site_list = self._pore.get_sites()
        self._site_in = [site_key for site_key, site_val in site_list.items() if site_val["type"]=="in"]
        self._site_ex = [site_key for site_key, site_val in site_list.items() if site_val["type"]=="ex"]

        # Save original surface si positions
        self._si_pos_in = [self._block.pos(site_key) for site_key, site_val in site_list.items() if site_val["type"]=="in"]
        self._si_pos_ex = [self._block.pos(site_key) for site_key, site_val in site_list.items() if site_val["type"]=="ex"]

        # Siloxane bridges
        if hydro[0]:
            self._siloxane(hydro[0], "in")
            self._site_in = [site_key for site_key, site_val in site_list.items() if site_val["type"]=="in"]
        if hydro[1] and res:
            self._siloxane(hydro[1], "ex")
            self._site_ex = [site_key for site_key, site_val in site_list.items() if site_val["type"]=="ex"]

        # Objectify grid
        non_grid = matrix.bound(1)+list(site_list.keys())
        bonded = matrix.bound(0, "gt")
        grid_atoms = [atom for atom in bonded if not atom in non_grid]
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

    def attach_special(self, mol, mount, axis, amount, scale=1, symmetry="point", is_proxi=True, is_rotate=False):
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
        is_proxi : bool, optional
            True to fill binding sites in proximity of filled binding site
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

            x = self._centroid[0]+coeff*self.height()/2
            y = self._centroid[1]
            z = start+dist*i

            pos_list.append([x, y, z])

        # Run attachment
        mols = self._pore.attach(mol, mount, axis, self._site_in, len(pos_list), self._normal_in, scale, pos_list=pos_list, is_proxi=is_proxi, is_random=False, is_rotate=is_rotate)

        # Add to sorting list
        for mol in mols:
            if not mol.get_short() in self._sort_list:
                self._sort_list.append(mol.get_short())


    ############
    # Analysis #
    ############
    def height(self):
        """Calculate true slit pore size after drilling and preperation. This
        is done by determining the mean value :math:`\\bar r` of the silicon
        distances :math:`r_i` of silicon :math:`i` towards the pore center

        .. math::

            \\bar r=\\frac1n\\sum_{i=1}^nr_i

        with the number of silicon atoms :math:`n`. The size is then

        .. math::

            h=2\\bar r=\\frac2n\\sum_{i=1}^nr_i.

        Returns
        -------
        height : float
            Pore size after preperation
        """
        # Calculate distance towards central axis of binding site silicon atoms
        r = [pms.geom.length(pms.geom.vector([pos[0], self._centroid[1], pos[2]], pos)) for pos in self._si_pos_in]

        # Calculate mean
        r_bar = sum(r)/len(r) if len(r)>0 else 0

        # Calculate diameter
        return 2*r_bar

    def roughness(self):
        """Calculate surface roughness. In the case of a slit pore the roughness
        is then determined by calculating the standard deviation of the binding
        site silicon atoms peaks and valleys.

        It is therefore enough to calculate the distances towards a specific
        surface, which in this case will be the central surface. The mean value
        :math:`\\bar r` of the silicon distances :math:`r_i` of silicon
        :math:`i` towards the pore center, is calculated by

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
        # Interior
        ## Calculate distance towards central axis of binding site silicon atoms
        r_in = [pms.geom.length(pms.geom.vector([pos[0], self._centroid[1], pos[2]], pos)) for pos in self._si_pos_in]

        # Exterior
        if self._res:
            ## Create molecules with exterior positions
            temp_mol = pms.Molecule()
            for pos in self._si_pos_ex:
                temp_mol.add("Si", pos)
            temp_mol.zero()
            size = temp_mol.get_box()[2]

        ## Calculate distance to boundary
        r_ex = [pos[2] if pos[2] < size/2 else abs(pos[2]-size) for pos in self._si_pos_ex]

        # Calculate mean
        r_bar_in = sum(r_in)/len(r_in) if len(r_in)>0 else 0
        r_bar_ex = sum(r_ex)/len(r_ex) if len(r_ex)>0 else 0

        # Calculate roughness
        r_q_in =  math.sqrt(sum([(r_i-r_bar_in)**2 for r_i in r_in])/len(r_in)) if len(r_in)>0 else 0
        r_q_ex =  math.sqrt(sum([(r_i-r_bar_ex)**2 for r_i in r_ex])/len(r_ex)) if len(r_ex)>0 else 0

        # Calculate square root roughness
        return {"in": r_q_in, "ex": r_q_ex}

    def volume(self):
        """Calculate pore volume. This is done by defining a new shape object
        with system sizes after pore preparation and using the volume function
        :func:`porems.shape.Cuboid.volume`.

        Returns
        -------
        volume : float
            Pore volume
        """
        return pms.Cuboid({"centroid": self._centroid, "central": [0, 0, 1], "length": self._box[2], "width": self._box[0], "height": self.height()}).volume()

    def surface(self):
        """Calculate pore surface and exterior surface. This is can be simply
        calculated by

        .. math::

            A = 2\\cdot x\\cdot z

        with block width :math:`x` and depth :math:`z`.

        Returns
        -------
        surface : dictionary
            Pore surface of interior and exterior
        """
        return {"in": 2*self._box[2]*self._box[0], "ex": self._box[0]*(self._box[1]-self.height())}

    def table(self, decimals=3):
        """Create properties as pandas table for easy viewing.

        Parameters
        ----------
        decimals : integer, optional
            Number of decimals to be rounded to

        Returns
        -------
        tables : DataFrame
            Pandas table of all properties
        """
        # Initialize
        form = "%."+str(decimals)+"f"

        # Define pore specific properties
        props = {}
        props["Pore height (nm)"] = [form%self.height(), " "]

        return self._table_base(props, decimals)


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
    hydro: list, optional
        Hydroxilation degree for interior and exterior of the pore in
        :math:`\\frac{\\mu\\text{mol}}{\\text{m}^2}`

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
    def __init__(self, size, diam, sep, res=5, hydro=[0, 0]):
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
        self._len_cyl = (self._pattern.get_size()[2]-sep-diam)/2

        # Create block
        self._block = self._pattern.get_block()
        self._block.set_name("pore")
        self._box = self._block.get_box()

        # Dice up block
        dice = pms.Dice(self._block, 0.4, True)
        matrix = pms.Matrix(dice.find_parallel(None, ["Si", "O"], [0.155-1e-2, 0.155+1e-2]))

        # Create pore object
        self._pore = pms.Pore(self._block, matrix)
        self._pore.set_name("pore")

        # Create exterior
        if res:
            self._pore.exterior(self._pattern.get_gap())

        # Carve out shape
        central = [0, 0, 1]

        self._centroid = {}
        self._centroid["block"] = self._block.centroid()
        self._centroid["cyl_l"] = self._centroid["block"][:2]+[0]
        self._centroid["cyl_r"] = self._centroid["block"][:2]+[self._pattern.get_size()[2]]
        self._centroid["sph_l"] = self._centroid["block"][:2]+[self._len_cyl]
        self._centroid["sph_r"] = self._centroid["block"][:2]+[self._pattern.get_size()[2]-self._len_cyl]

        self._shape = {}
        self._shape["cyl_l"] = pms.Cylinder({"centroid": self._centroid["cyl_l"], "central": central, "length": self._len_cyl*2, "diameter": diam-0.5})
        self._shape["cyl_r"] = pms.Cylinder({"centroid": self._centroid["cyl_r"], "central": central, "length": self._len_cyl*2, "diameter": diam-0.5})
        self._shape["sph_l"] = pms.Sphere({"centroid": self._centroid["sph_l"], "central": central, "diameter": diam-0.5})
        self._shape["sph_r"] = pms.Sphere({"centroid": self._centroid["sph_r"], "central": central, "diameter": diam-0.5})

        del_list = []
        for shape in self._shape.values():
            del_list.extend([atom_id for atom_id, atom in enumerate(self._block.get_atom_list()) if shape.is_in(atom.get_pos())])
        matrix.strip(del_list)

        self._pore_shape = "CAPSULE"

        # Prepare pore surface
        self._pore.prepare()

        # Determine sites
        self._pore.sites()
        site_list = self._pore.get_sites()
        self._site_in = [site_key for site_key, site_val in site_list.items() if site_val["type"]=="in"]
        self._site_ex = [site_key for site_key, site_val in site_list.items() if site_val["type"]=="ex"]

        # Save original surface si positions
        self._si_pos_in = [self._block.pos(site_key) for site_key, site_val in site_list.items() if site_val["type"]=="in"]
        self._si_pos_ex = [self._block.pos(site_key) for site_key, site_val in site_list.items() if site_val["type"]=="ex"]

        # Siloxane bridges
        if hydro[0]:
            self._siloxane(hydro[0], "in")
            self._site_in = [site_key for site_key, site_val in site_list.items() if site_val["type"]=="in"]
        if hydro[1]:
            self._siloxane(hydro[1], "ex")
            self._site_ex = [site_key for site_key, site_val in site_list.items() if site_val["type"]=="ex"]

        # Objectify grid
        non_grid = matrix.bound(1)+list(site_list.keys())
        bonded = matrix.bound(0, "gt")
        grid_atoms = [atom for atom in bonded if not atom in non_grid]
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


    ############
    # Analysis #
    ############
    def diameter(self):
        """Calculate true cylinder diameter after drilling and preperation. This
        is done by determining the mean value :math:`\\bar r` of the silicon
        distances :math:`r_i` of silicon :math:`i` towards the pore center

        .. math::

            \\bar r=\\frac1n\\sum_{i=1}^nr_i

        with the number of silicon atoms :math:`n`. Hereby the silicon atoms are
        only of the cylidric parts of the capsules. The diameter is then

        .. math::

            d=2\\bar r=\\frac2n\\sum_{i=1}^nr_i.

        Returns
        -------
        diameter : float
            Pore diameter after preperation
        """
        # Determine sites in the cylindric part of the pore
        pos_list = []
        for pos in self._si_pos_in:
            if pos[2] <= self._len_cyl:
                pos_list.append(pos)
            elif pos[2] >= self._pattern.get_size()[2]-self._len_cyl:
                pos_list.append(pos)

        # Calculate distance towards central axis of binding site silicon atoms
        r = [pms.geom.length(pms.geom.vector([self._centroid["block"][0], self._centroid["block"][1], pos[2]], pos)) for pos in pos_list]

        # Calculate mean
        r_bar = sum(r)/len(r) if len(r)>0 else 0

        # Calculate diameter
        return 2*r_bar

    def roughness(self):
        """Calculate surface roughness. In the case of a capsule pore one can
        visualize pulling the cylindric section pore apart, thus flattening the
        interior surface. The roughness is then determined by calculating the
        standard deviation of the binding site silicon atoms peaks and valleys.

        It is therefore enough to calculate the distances towards a specific
        axis, which in this case will be the central axis. The mean value
        :math:`\\bar r` of the silicon distances :math:`r_i` of silicon
        :math:`i` towards the pore center, is calculated by

        .. math::

            \\bar r=\\frac1n\\sum_{i=1}^nr_i

        with the number of silicon atoms :math:`n`. Hereby only the silicon
        atoms within the cyldirc part are considered. This mean value is used in
        the square root roughness calculation

        .. math::

            R_q = \\sqrt{\\frac1n\\sum_{i=1}^n\\|r_i-\\bar r\\|^2}.

        Returns
        -------
        roughness : float
            Surface roughness
        """
        # Interior
        ## Determine sites in the cylindric part of the pore
        pos_list = []
        for pos in self._si_pos_in:
            if pos[2] <= self._len_cyl:
                pos_list.append(pos)
            elif pos[2] >= self._pattern.get_size()[2]-self._len_cyl:
                pos_list.append(pos)

        ## Calculate distance towards central axis of binding site silicon atoms
        r_in = [pms.geom.length(pms.geom.vector([self._centroid["block"][0], self._centroid["block"][1], pos[2]], pos)) for pos in pos_list]

        # Exterior
        ## Create molecules with exterior positions
        if self._res:
            temp_mol = pms.Molecule()
            for pos in self._si_pos_ex:
                temp_mol.add("Si", pos)
            temp_mol.zero()
            size = temp_mol.get_box()[2]

        ## Calculate distance to boundary
        r_ex = [pos[2] if pos[2] < size/2 else abs(pos[2]-size) for pos in self._si_pos_ex]

        # Calculate mean
        r_bar_in = sum(r_in)/len(r_in) if len(r_in)>0 else 0
        r_bar_ex = sum(r_ex)/len(r_ex) if len(r_ex)>0 else 0

        # Calculate roughness
        r_q_in =  math.sqrt(sum([(r_i-r_bar_in)**2 for r_i in r_in])/len(r_in)) if len(r_in)>0 else 0
        r_q_ex =  math.sqrt(sum([(r_i-r_bar_ex)**2 for r_i in r_ex])/len(r_ex)) if len(r_ex)>0 else 0

        # Calculate square root roughness
        return {"in": r_q_in, "ex": r_q_ex}

    def volume(self):
        """Calculate pore volume. This is done by defining a new shape object
        with system sizes after pore preparation and using the volume function
        :func:`porems.shape.Cylinder.volume`.

        Returns
        -------
        volume : float
            Pore volume
        """
        diam = self.diameter()

        cylinder = pms.Cylinder({"centroid": [0, 0, 0], "central": [0, 0, 1], "length": self._len_cyl*2, "diameter": diam})
        sphere = self._shape["sph_l"] = pms.Sphere({"centroid": [0, 0, 0], "central": [0, 0, 1], "diameter": diam})

        return cylinder.volume() + sphere.volume()

    def surface(self):
        """Calculate pore surface and exterior surface. This is done by defining
        a new shape object with system sizes after pore preparation and using
        the surface function :func:`porems.shape.Cylinder.surface`.

        Returns
        -------
        surface : dictionary
            Pore surface of interior and exterior
        """
        diam = self.diameter()

        cylinder = pms.Cylinder({"centroid": [0, 0, 0], "central": [0, 0, 1], "length": self._len_cyl*2, "diameter": diam})
        sphere = self._shape["sph_l"] = pms.Sphere({"centroid": [0, 0, 0], "central": [0, 0, 1], "diameter": diam})

        surf_in = cylinder.surface() + sphere.surface()
        surf_ex = 2*(self._box[0]*self._box[1]-math.pi*(diam/2)**2)

        return {"in": surf_in, "ex": surf_ex}

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
        form = "%."+str(decimals)+"f"

        # Define pore specific properties
        props = {}
        props["Pore diameter (nm)"] = [form%self.diameter(), " "]
        props["Cavity separation distance(nm)"] = [form%self._sep, " "]

        return self._table_base(props, decimals)


class PoreAmorphCylinder(PoreSystem):
    """This class carves a cylindric pore system out of an amorph
    :math:`\\beta`-cristobalite block from
    `Vink et al. <http://doi.org/10.1103/PhysRevB.67.245201>`_ with dimensions
    [9.605, 9.605, 9.605] (x, y, z).

    Parameters
    ----------
    diam : float
        Cylinder diameter
    res : float, optional
        Reservoir size on each side
    hydro: list, optional
        Hydroxilation degree for interior and exterior of the pore in
        :math:`\\frac{\\mu\\text{mol}}{\\text{m}^2}`

    Examples
    --------
    Following example generates a cylindrical pore with a diameter of 4nm,
    reservoirs of 5nm on each side and a surface functionalized with TMS

    .. code-block:: python

        import porems as pms

        pore = pms.PoreAmorphCylinder(4, 5)

        pore.attach(pms.gen.tms(), 0, [0, 1], 100, "in")
        pore.attach(pms.gen.tms(), 0, [0, 1], 20, "ex")

        pore.finalize()

        pore.store("output/")
    """
    def __init__(self, diam, res=5, hydro=[0, 0]):
        # Call super class
        super(PoreAmorphCylinder, self).__init__()

        # Initialize
        self._diam = diam
        self._res = res
        self._sort_list = ["OM", "SI"]

        # Build pattern
        pattern = pms.BetaCristobalit()
        pattern.generate([2, 2, 2], "z")
        pattern._structure = pms.Molecule(inp=os.path.split(__file__)[0]+"/templates/amorph.gro")
        pattern._size = [2.014, 1.751, 2.468]

        # Create block
        self._block = pattern.get_block()
        self._block.set_name("pore")
        self._box = self._block.get_box()

        # Dice up block
        dice = pms.Dice(self._block, 0.4, True)
        matrix = pms.Matrix(dice.find_parallel(None, ["Si", "O"], [0.160-0.02, 0.160+0.02]))
        matrix.split(57790, 2524)

        # Create pore object
        self._pore = pms.Pore(self._block, matrix)
        self._pore.set_name("pore")

        # Create exterior
        if res:
            self._pore.exterior(pattern.get_gap())

        # Carve out shape
        self._centroid = self._block.centroid()
        central = [0, 0, 1]
        self._cylinder = pms.Cylinder({"centroid": self._centroid, "central": central, "length": 9.605, "diameter": diam-0.5})  # Preperation precaution
        del_list = [atom_id for atom_id, atom in enumerate(self._block.get_atom_list()) if self._cylinder.is_in(atom.get_pos())]
        matrix.strip(del_list)
        self._pore_shape = "CYLINDER"

        # Prepare pore surface
        self._pore.prepare()

        # Determine sites
        self._pore.sites()
        site_list = self._pore.get_sites()
        self._site_in = [site_key for site_key, site_val in site_list.items() if site_val["type"]=="in"]
        self._site_ex = [site_key for site_key, site_val in site_list.items() if site_val["type"]=="ex"]

        # Save original surface si positions
        self._si_pos_in = [self._block.pos(site_key) for site_key, site_val in site_list.items() if site_val["type"]=="in"]
        self._si_pos_ex = [self._block.pos(site_key) for site_key, site_val in site_list.items() if site_val["type"]=="ex"]

        # Siloxane bridges
        if hydro[0]:
            self._siloxane(hydro[0], "in", [0, 0.51])
            self._site_in = [site_key for site_key, site_val in site_list.items() if site_val["type"]=="in"]
        if hydro[1]:
            self._siloxane(hydro[1], "ex", [0, 0.51])
            self._site_ex = [site_key for site_key, site_val in site_list.items() if site_val["type"]=="ex"]

        # Objectify grid
        non_grid = matrix.bound(1)+list(site_list.keys())
        bonded = matrix.bound(0, "gt")
        grid_atoms = [atom for atom in bonded if not atom in non_grid]
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

    def attach_special(self, mol, mount, axis, amount, scale=1, symmetry="point", is_proxi=True, is_rotate=False):
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
        is_proxi : bool, optional
            True to fill binding sites in proximity of filled binding site
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
        mols = self._pore.attach(mol, mount, axis, self._site_in, len(pos_list), self._normal_in, scale, pos_list=pos_list, is_proxi=is_proxi, is_random=False, is_rotate=is_rotate)

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
        r = [pms.geom.length(pms.geom.vector([self._centroid[0], self._centroid[1], pos[2]], pos)) for pos in self._si_pos_in]

        # Calculate mean
        r_bar = sum(r)/len(r) if len(r)>0 else 0

        # Calculate diameter
        return 2*r_bar

    def roughness(self):
        """Calculate surface roughness. In the case of a cylindric pore one can
        visualize pulling the pore apart, thus flattening the interior surface.
        The roughness is then determined by calculating the standard deviation
        of the binding site silicon atoms peaks and valleys.

        It is therefore enough to calculate the distances towards a specific
        axis, which in this case will be the central axis. The mean value
        :math:`\\bar r` of the silicon distances :math:`r_i` of silicon
        :math:`i` towards the pore center, is calculated by

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
        # Interior
        ## Calculate distance towards central axis of binding site silicon atoms
        r_in = [pms.geom.length(pms.geom.vector([self._centroid[0], self._centroid[1], pos[2]], pos)) for pos in self._si_pos_in]

        # Exterior
        if self._res:
            ## Create molecules with exterior positions
            temp_mol = pms.Molecule()
            for pos in self._si_pos_ex:
                temp_mol.add("Si", pos)
            temp_mol.zero()
            size = temp_mol.get_box()[2]

        ## Calculate distance to boundary
        r_ex = [pos[2] if pos[2] < size/2 else abs(pos[2]-size) for pos in self._si_pos_ex]

        # Calculate mean
        r_bar_in = sum(r_in)/len(r_in) if len(r_in)>0 else 0
        r_bar_ex = sum(r_ex)/len(r_ex) if len(r_ex)>0 else 0

        # Calculate roughness
        r_q_in =  math.sqrt(sum([(r_i-r_bar_in)**2 for r_i in r_in])/len(r_in)) if len(r_in)>0 else 0
        r_q_ex =  math.sqrt(sum([(r_i-r_bar_ex)**2 for r_i in r_ex])/len(r_ex)) if len(r_ex)>0 else 0

        # Calculate square root roughness
        return {"in": r_q_in, "ex": r_q_ex}

    def volume(self):
        """Calculate pore volume. This is done by defining a new shape object
        with system sizes after pore preparation and using the volume function
        :func:`porems.shape.Cylinder.volume`.

        Returns
        -------
        volume : float
            Pore volume
        """
        return pms.Cylinder({"centroid": self._centroid, "central": [0, 0, 1], "length": self._box[2], "diameter": self.diameter()}).volume()

    def surface(self):
        """Calculate pore surface and exterior surface. This is done by defining
        a new shape object with system sizes after pore preparation and using
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

    def table(self, decimals=3):
        """Create properties as pandas table for easy viewing.

        Parameters
        ----------
        decimals : integer, optional
            Number of decimals to be rounded to

        Returns
        -------
        tables : DataFrame
            Pandas table of all properties
        """
        # Initialize
        form = "%."+str(decimals)+"f"

        # Define pore specific properties
        props = {}
        props["Pore diameter (nm)"] = [form%self.diameter(), " "]

        return self._table_base(props, decimals)
