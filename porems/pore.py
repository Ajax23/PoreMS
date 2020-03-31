################################################################################
# Pore Class                                                                   #
#                                                                              #
"""Function for editing a pore surface."""
################################################################################


from collections import Counter


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

        self._block = block
        self._matrix = matrix


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

    def sites(self, exterior=None):
        """Create binding site dictionary of the format

        .. math::

            \\boldsymbol{B}=\\begin{Bmatrix}
                si_1:&
                \\begin{Bmatrix}
                    \\text{"o"}: \\begin{pmatrix}o_{1,1},\\dots o_{1,n}\\end{pmatrix}&
                    \\text{"type"}: \\text{"in"/"ex"}&
                    \\text{"state"}: 0/1/2
                \\end{Bmatrix}\\\\
                \\vdots
            \\end{bmatrix}

        with entries

        * **o** - Unsaturated oxygen atoms bound to silicon atom
        * **type** - Exterior "ex" or interior "in" binding site
        * **state** - available - 0, used - 1, unused but deactivated - 2

        Parameters
        ----------
        exterior : list, None, optional
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
            if exterior is not None:
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
            data["state"] = 0


    #######################
    # Molecule Attachment #
    #######################



    ###############
    # Final Edits #
    ###############



    ###########
    # Analyze #
    ###########



    ##################
    # Setter Methods #
    ##################



    ##################
    # Getter Methods #
    ##################
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
