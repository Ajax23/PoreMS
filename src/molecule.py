################################################################################
# Molecule Class                                                               #
#                                                                              #
"""All necessary function for creating and editing molecules."""
################################################################################


import math
import config

from build.database import Database


class Molecule:
    """This class defines a molecule object, which is basically a data matrix
    that contains atom information of positions and type.

    Be Molecule :math:`\\boldsymbol{D}\\in\\mathbb{R}^{(n+1)\\times m}` with spatial dimensions
    :math:`i=1,\\dots,n` and atom number :math:`j=1,\\dots,m`.
    The data matrix then has following structure

    .. math::

        \\boldsymbol{D}=\\begin{bmatrix}
        \\begin{pmatrix}d_{11}\\\\d_{12}\\\\\\vdots\\\\d_{1m}\\end{pmatrix}&
        \\begin{pmatrix}d_{21}\\\\d_{22}\\\\\\vdots\\\\d_{2m}\\end{pmatrix}&
        \\dots&
        \\begin{pmatrix}d_{n1}\\\\d_{n2}\\\\\\vdots\\\\d_{nm}\\end{pmatrix}&
        \\begin{pmatrix}d_{t1}\\\\d_{t2}\\\\\\vdots\\\\d_{tm}\\end{pmatrix}
        \\end{bmatrix}

    with atom data :math:`d_{ij}` of atom :math:`j` and chemical atom type entry :math:`t:=n+1`.

    The requested dimension is defined in the parameter module.
    Most of the following functions have been written for :math:`n`-dimensions.
    However they only have been tested for the three-dimnesional case.

    Functions have been provided for editing, moving and transforming the
    entries of the data matrix either for the whole matrix or for specified atoms.
    Private functions are solely needed for internal puposes,
    whereas the public functions are inteded to be used for the editingself.

    For example molecule consider class extensions **Catalysator** and **Educt**.

    Parameters
    ----------
    name : str
        Molecule name
    short : str
        Molecule short name
    inp : None, str, list
        None for empty Molecule, string to read a molecule from a
        specified filelink or a list of either molecules to concatenate these
        into one object, or a data list in the earlier discussed format
    """
    def __init__(self,inp=None,name="Molecule",short="MOL"):
        # Initialize
        self.__box    = None
        self.__masses = None
        self.__mass   = None
        self.__com    = None
        self.__write  = [self]

        self.__dim    = 3
        self.__charge = 0

        self.__name   = name
        self.__short  = short

        # Check data input
        if inp is None:
            self.__data = [[] for i in range(self.__dim+1)]
        else:
            # Read from file
            if isinstance(inp,str):
                self.__data = self.__read(inp,inp.split(".")[-1].upper())
            # Concat multiple molecules
            elif isinstance(inp,list):
                # Data list is provided
                if(isinstance(inp[0],list)): self.__data = inp
                # List of molecules is provided
                else: self.__data = self.__concat(inp)


    ################################
    # Private Methods - Management #
    ################################
    def __read(self,file_name,file_type):
        """Read a molecule from a file. Currently only **GRO**, **PFB** and
        **MOL2** files are supported.

        Parameters
        ----------
        file_name : str
            Link to requested file
        file_type : str
            Int for types list id or str for file extension name

        Returns
        -------
        data : list
            data matrix
        """
        # Process input
        if not file_type in ["GRO","PDB","MOL2"]:
            print("Unsupported filetype.")
            return

        # Read molecules
        data = []
        with open(file_name,"r") as file_in:
            line_idx = 0
            for line in file_in:
                line_val = line.split()

                # Gro file
                if file_type=="GRO":
                    if line_idx>0 and len(line_val)>3:
                        coord = [float(line_val[i]) for i in range(3,5+1)]
                        coord.append(''.join([i for i in line_val[1] if not i.isdigit()]))
                        data.append(coord)

                # Pdb file
                elif file_type=="PDB":
                    if line_val[0] in ["ATOM","HETATM"]:
                        coord = [float(line_val[i])/10 for i in range(6,7+1)]
                        coord.append(line_val[11])
                        data.append(coord)

                # Mol2 file
                elif file_type=="MOL2":
                    if len(line_val)>8:
                        coord = [float(line_val[i])/10 for i in range(2,4+1)]
                        coord.append(''.join([i for i in line_val[1] if not i.isdigit()]))
                        data.append(coord)

                # Running variable
                line_idx += 1

        # Transform to column
        return utility.column(data)

    def __concat(self,mol_list):
        """Concatenate a molecule list into one molecule object.

        Parameters
        ----------
        mol_list : list
            List of molecule objects to be concatenated

        Returns
        -------
        data : list
            data matrix
        """
        # Initialize
        data = [[] for x in range(self.__dim+1)]

        # Concatenate
        for mol in mol_list:
            # Molecules
            mol_data = mol.getData()
            for i in range(len(mol_data)):
                for j in range(len(mol_data[i])):
                    data[i].append(mol_data[i][j])

        return data

    def __temp(self,atoms):
        """Create a temporary molecule of specified atom ids.

        Parameters
        ----------
        atoms : list
            List of atoms to be duplicated

        Returns
        -------
        mol : Molecule
            Molecule object
        """
        # Create temp matrix
        temp = [[] for x in range(self.__dim+1)]

        # Fill matrix
        for i in range(len(atoms)):
            for j in range(self.__dim+1):
                temp[j].append(self.__data[j][atoms[i]])

        # Create molecule
        return Molecule(inp=temp)

    def __update(self,data,atoms):
        """Exchange data of specified atom with given data.

        Parameters
        ----------
        data : list
            New data list
        atoms : list
            List of atoms to be updated
        """
        for i in range(len(atoms)):
            for j in range(self.__dim):
                self.__data[j][atoms[i]] = data[j][i]

    def __append(self,mol):
        """Append a given molecule to the current object.

        Parameters
        ----------
        mol : Molecule
            Molecule object
        """
        temp = mol.getData()
        for i in range(len(temp)):
            for j in range(len(temp[i])):
                self.__data[i].append(temp[i][j])


    ##############################
    # Private Methods - Geometry #
    ##############################
    def __dotproduct(self,vec1,vec2):
        """Calculate the dot-product of two vectors
        :math:`\\boldsymbol{a},\\boldsymbol{b}\\in\\mathbb{R}^n`

        .. math::

            dot(\\boldsymbol{a},\\boldsymbol{b})=
            \\begin{pmatrix}a_1\\\\\\vdots\\\\a_n\\end{pmatrix}\\cdot
            \\begin{pmatrix}b_1\\\\\\vdots\\\\b_n\\end{pmatrix}=
            a_1\\cdot b_1+a_2\\cdot b_2+\\dots+a_n\\cdot b_n.

        Parameters
        ----------
        vec1 : list
            First vector
        vec2 : list
            Second vector

        Returns
        -------
        dot : float
            Dotproduct value
        """
        return sum((a*b) for a, b in zip(vec1,vec2))

    def __length(self,vec):
        """Calculate the length of a vector :math:`\\boldsymbol{a}\\in\\mathbb{R}^n`

        .. math::

            length(\\boldsymbol{a})=|\\boldsymbol{a}|=\\sqrt{\\boldsymbol{a}\cdot\\boldsymbol{a}}

        Parameters
        ----------
        vec : list
            Input vector

        Returns
        -------
        length : float
            Vector length
        """
        return math.sqrt(self.__dotproduct(vec,vec))

    def __vector(self,inpA,inpB):
        """Calculate the vector beween to two positions
        :math:`\\boldsymbol{a},\\boldsymbol{b}\\in\\mathbb{R}^n`

        .. math::

            vec(\\boldsymbol{a},\\boldsymbol{b})=\\begin{pmatrix}b_1-a_1\\\\\\vdots\\\\b_n-a_n\\end{pmatrix}

        The two inputs can either be atom indices or to vectorial positions.

        Parameters
        ----------
        inpA : int, list
            Either an atom id or a position vector
        inpB : int, list
            Either an atom id or a position vector

        Returns
        -------
        vector : list
            Bond vector
        """
        # Initialize
        dim = self.__dim
        vec = []

        # Process input
        if isinstance(inpA,int) and isinstance(inpB,int):
            inpA = self.pos(inpA)
            inpB = self.pos(inpB)
        elif not isinstance(inpA,list) and isinstance(inpB,list):
            print("Vector: Wrong input...")
            return

        # Check dimensions
        if not len(inpA)==len(inpB) or not len(inpA)==dim:
            print("Vector: Wrong dimensions...")
            return

        # Calculate vector
        for i in range(dim):
            vec.append(inpB[i]-inpA[i])
        return vec

    def __unit(self,vec):
        """Transform a vector :math:`\\boldsymbol{a}\\in\\mathbb{R}^n` into a unit vactor

        .. math::

            unit(\\boldsymbol{a})=\\frac{\\boldsymbol{a}}{|\\boldsymbol{a}|}

        Parameters
        ----------
        pos1 : list
            First vector
        pos1 : list
            Second vector

        Returns
        -------
        vec : list
            Vector
        """
        vec    = vec[:]
        length = self.__length(vec)
        for i in range(self.__dim):
            vec[i] = vec[i]/length if not length==0 else vec[i]

        return vec

    def __cross(self,vecA,vecB):
        """Calculate the crossproduct of two three-dimensional vectors
        :math:`\\boldsymbol{a},\\boldsymbol{b}\\in\\mathbb{R}^3`

        .. math::

            cross(\\boldsymbol{a},\\boldsymbol{b})=\\begin{pmatrix}
            a_2\\cdot b_3-a_3\\cdot b_2\\\\
            a_3\\cdot b_1-a_1\\cdot b_4\\\\
            a_1\\cdot b_2-a_2\\cdot b_1
            \\end{pmatrix}

        Parameters
        ----------
        vecA : list
            Vector a
        vecB : list
            Vector b

        Returns
        -------
        vec : list
            Crossproduct Vector
        """
        vec = []
        vec.append(vecA[1]*vecB[2]-vecA[2]*vecB[1])
        vec.append(vecA[2]*vecB[0]-vecA[0]*vecB[2])
        vec.append(vecA[0]*vecB[1]-vecA[1]*vecB[0])

        return vec

    def __anglePolar(self,pos,isDeg=False):
        """Calculate the polar angle of a position vector
        :math:`\\boldsymbol{a}\\in\\mathbb{R}^3`, which is the angle of the
        x-axis towards the reflected position vector on the x-y-plane

        .. math::

            polar(\\boldsymbol{a})=\\arctan2(x,y)\\left\\{
            \\begin{array}{ll}
            \\tan^{-1}\\left(\\frac{y}{x}\\right)&x>0\\\\
            \\tan^{-1}\\left(\\frac{y}{x}\\right)+\\pi&x<0,y>0\\\\
            \\pm\\pi&x<0,y=0\\\\
            \\tan^{-1}\\left(\\frac{y}{x}\\right)-\\pi&x<0,y<0\\\\
            +\\frac{\\pi}{2}&x=0,y>0\\\\
            -\\frac{\\pi}{2}&x=0,y<0
            \\end{array}
            \\right.

        with :math:`x` as the first vector entry and :math:`y` as the second.

        Parameters
        ----------
        pos : list
            Position vector a
        isDeg : bool
            True if the output should be in degree

        Returns
        -------
        angle : float
            Polar angle
        """
        try:
            angle = math.atan2(pos[1],pos[0])
        except:
            angle = math.atan(0)

        return angle*180/math.pi if isDeg else angle

    def __angleAzi(self,pos,isDeg=False):
        """Calculate the azimuthal angle of a position vector
        :math:`\\boldsymbol{a}\\in\\mathbb{R}^3`, which is the angle of the
        position vector towards the x-y-plane

        .. math::

            azimut(\\boldsymbol{a})=\\cos^{-1}\\frac{y}{|\\boldsymbol{a}|}

        with :math:`y` as the second vector entry.

        Parameters
        ----------
        pos : list
            Position vector a
        isDeg : bool
            True if the output should be in degree

        Returns
        -------
        angle : float
            Azimuthal angle
        """
        try:
            angle = math.acos(pos[2]/self.__length(pos))
        except:
            angle = math.acos(0)

        return angle*180/math.pi if isDeg else angle

    def __axis(self,inp):
        """Return the unit-vector of the main axes. Input is either integer or string with

        * 1 or "x" for the x-axis
        * 2 or "y" for the y-axis
        * 3 or "z" for the z-axis

        Parameters
        ----------
        inp : int, str
            Axis type input

        Returns
        -------
        vec : list
            Unit vector
        """
        # Error message
        axisError = "Wrong axis definition..."

        # Process input
        if   isinstance(inp,str):
            if   inp=="x":
                axis = 1
            elif inp=="y":
                axis = 2
            elif inp=="z":
                axis = 3
            else:
                return axisError
        elif isinstance(inp,int):
            if inp==1 or inp==2 or inp==3:
                axis = inp
            else:
                return axisError
        else:
            return axisError

        # Return vector
        return [1 if i==axis-1 else 0 for i in range(self.__dim)]

    def __focal(self):
        """Calculate the positional focal-point

        .. math::

            focal=\\begin{pmatrix}f_1&f_2&\\dots&f_n\\end{pmatrix}^T

        with

        .. math::

            f_i=\\frac{1}{m}\\sum_j^m d_{ij}

        Returns
        -------
        focal : list
            Focal point
        """
        # Initialize
        data = self.__data

        # Calculate focal point
        return [sum(data[i])/len(data[i]) for i in range(self.__dim)]

    def __angle(self,v1,v2,isDeg=True):
        """Calculate the angle between two vectors
        :math:`\\boldsymbol{a}\\in\\mathbb{R}^n`

        .. math::

            angle=\\cos^{-1}\\frac{\\boldsymbol{a}\cdot\\boldsymbol{b}}
            {|\\boldsymbol{a}||\\boldsymbol{a}|}

        Parameters
        ----------
        v1 : list
            First vector
        v2 : list
            Second vector
        isDeg : bool
            True if the output should be in degree

        Returns
        -------
        angle : float
            Angle
        """
        angle = math.acos(self.__dotproduct(v1,v2)/(self.__length(v1)*self.__length(v2)))

        return angle*180/math.pi if isDeg else angle

    def __rotate(self,data,axis,angle,isDeg):
        """Rotate a vector :math:`\\boldsymbol{a}\\in\\mathbb{R}^3`
        along an axis :math:`\\boldsymbol{b}\\in\\mathbb{R}^3` with angle
        :math:`\\alpha\\in\\mathbb{R}`.
        Input for the axis is either a vector or the input for
        function :func:`_vec`.
        The rotation is performed using the rotation-matrix

        .. math::

            \\boldsymbol{R}_\\boldsymbol{n}(\\alpha)=\\begin{pmatrix}
            n_1^2 (1-\\cos\\alpha)+   \\cos\\alpha&n_1n_2(1-\\cos\\alpha)-n_3\\sin\\alpha&n_1n_3(1-\\cos\\alpha)+n_2\\sin\\alpha\\\\
            n_2n_1(1-\\cos\\alpha)+n_3\\sin\\alpha&n_2^2 (1-\\cos\\alpha)+   \\cos\\alpha&n_2n_3(1-\\cos\\alpha)-n_1\\sin\\alpha\\\\
            n_3n_1(1-\\cos\\alpha)-n_2\\sin\\alpha&n_3n_2(1-\\cos\\alpha)+n_1\\sin\\alpha&n_3^2 (1-\\cos\\alpha)+   \\cos\\alpha
            \\end{pmatrix}


        where :math:`n_i` are the entries for the unit vector
        :math:`\\boldsymbol{n}` of the axis. The new coordinates :math:`\\boldsymbol{c}`
        are then calculated using a matrix-vector multiplication

        .. math::

            \\boldsymbol{c}=\\boldsymbol{R}_\\boldsymbol{n}\\boldsymbol{a}.

        Parameters
        ----------
        data : list
            Vector a
        axis : int, str, list
            Axis b
        angle : float
            Angle
        isDeg : bool
            True if the input is in degree

        Returns
        -------
        coord : list
            Vector c as the result of the rotation
        """
        # Initialize
        dim = self.__dim

        # Angle
        angle = angle*math.pi/180 if isDeg else angle

        # Set vector
        if isinstance(axis,list):
            if   len(axis)==dim:
                n = axis
            elif len(axis)==2:
                n = self.__vector(axis[0],axis[1])
            else:
                print("Wrong vector dimensions.")
                return
        else:
            n = self.__axis(axis)
            if isinstance(n,str):
                print(n)
                return

        # Unity
        n = self.__unit(n)

        # Define rotation matrix
        n1 = n[0]
        n2 = n[1]
        n3 = n[2]

        c  = math.cos(angle)
        s  = math.sin(angle)

        r  = [[n1*n1*(1.-c) +    c, n1*n2*(1.-c) - n3*s,  n1*n3*(1.-c) + n2*s],
              [n2*n1*(1.-c) + n3*s, n2*n2*(1.-c) +    c,  n2*n3*(1.-c) - n1*s],
              [n3*n1*(1.-c) - n2*s, n3*n2*(1.-c) + n1*s,  n3*n3*(1.-c) +    c]]

        # Rotate
        coord = []
        for i in range(dim):
            coord.append(data[0]*r[i][0]+data[1]*r[i][1]+data[2]*r[i][2])

        return coord

    def __boxSize(self):
        """Calculate the boxsize of the current molecule. This is done by
        determining the maximal coordinate value of all atoms in all dimensions

        .. math::

            \\boldsymbol{b}=\\begin{pmatrix}\\max(\\boldsymbol{d}_1)&max(\\boldsymbol{d}_1)&\\dots&max(\\boldsymbol{d}_n)\\end{pmatrix}^T

        where :math:`\\boldsymbol{d}_i` is the dimension-vector of the data matrix.

        Returns
        -------
        box : list
            Box length of the current molecule
        """
        data = self.__data
        dim  = self.__dim

        return [max(data[i]) if max(data[i])>0 else 0.001 for i in range(dim)]


    ##############################
    # Public Methods - Transform #
    ##############################
    def translate(self,vec):
        """Translate data matrix :math:`\\boldsymbol{D}` along a vector
        :math:`\\boldsymbol{a}\\in\\mathbb{R}^n`.

        .. math::

            \\boldsymbol{D}_{translate}=
            \\boldsymbol{D}+\\boldsymbol{a}=
            \\begin{pmatrix}
            \\boldsymbol{d}_1+a_1&\\boldsymbol{d}_2+a_2&\\dots&\\boldsymbol{d}_n+a_n&\\boldsymbol{d}_t
            \\end{pmatrix}

        Parameters
        ----------
        vec : list
            Vector a
        """
        # Initialize
        data = self.__data

        # Translate coordinates
        for i in range(self.__dim):
            for j in range(len(data[i])):
                data[i][j] += vec[i]

    def rotate(self,axis,angle,isDeg=True):
        """Rotate data matrix :math:`\\boldsymbol{D}` around an axis
        :math:`\\boldsymbol{a}\\in\\mathbb{R}^3` with angle
        :math:`\\alpha\\in\\mathbb{R}` using rotation function :func:`_rotate`.

        Parameters
        ----------
        axis : int, str, list
            Axis
        angle : float
            Angle
        isDeg : bool
            True if the input is in degree
        """
        # Initialize
        dim     = self.__dim
        data    = self.__data
        dataLen = len(data[0])

        # Rotate
        for i in range(dataLen):
            coord = [0 for i in range(dim)]
            for j in range(dim):
                coord[j] += data[j][i]

            coord = self.__rotate(coord,axis,angle,isDeg)
            for j in range(dim):
                data[j][i] = coord[j]

    def put(self,atom,pos):
        """Change the position of an atom to a given vector
        :math:`\\boldsymbol{a}\\in\\mathbb{R}^n`.

        Parameters
        ----------
        atom : int, str
            Atom-id whose position will be changed, use "last" for the newest atom
        pos : list
            New position vector
        """
        # User input
        if atom=="last":
            atom = len(self.__data[0])-1

        # Change position
        for i in range(len(pos)):
            self.__data[i][atom] = pos[i]

    def move(self,atom,pos):
        """Move whole the molecule to a new position, where the dragging point
        is a given atom that is moved to a specified position.

        Parameters
        ----------
        atom : int, str
            Main atom-id whose position will be changed, use "last" for the newest atom
        pos : list
            New position vector
        """
        # Calculate vector
        if atom=="last":
            atom = len(self.__data[0])-1

        vec = self.__vector(self.pos(atom),pos)

        # Translate
        self.translate(vec)

    def zero(self,pos=[0,0,0]):
        """Move whole the molecule, so that the minimal coordinate
        between all atoms is zero in all dimensions, or rather the values of the
        position variable *pos*. This function is basically setting
        the zero point of the coordinate system to *pos*.

        Parameters
        ----------
        pos : list
            Vector of the zero point of the coordinate system

        Returns
        -------
        vec : list
            Vector used for the translation
        """
        # Calculate translation vector
        vec = [pos[i]-min(self.__data[i]) for i in range(self.__dim)]

        # Resett box size
        self.__box = None

        # Translate molecule
        self.translate(vec)

        return vec

    def box(self,size,isExtend=True):
        """Create a molecule box around a centered molecule, where the
        coordinate systems zero point is set to the positional focal point of
        the molecule. The resulting box has the edge length of the given size
        value. If *isExtend* is set to True, then the maximal molecule
        coordinate of all dimensions is added to *size*.

        Parameters
        ----------
        size : float
            Box size
        isExtend : bool
            Should the molecule size be added to the box-size
        """
        # Initialize
        data = self.__data
        dim  = self.__dim

        # Get focal point
        focal = self.__focal()

        # Calculate box size
        boxFocal  = [f+size for f in focal]
        self.__box = [max(boxFocal) if isExtend else size for f in range(3)]

        # Move molecule to box center
        self.add   ("R",focal)
        self.move  ("last",[b/2 for b in self.__box])
        self.delete("last")


    #########################
    # Public Methods - Edit #
    #########################
    def partMove(self,bond,atoms,length,vec=None):
        """Change the length of a specified bond. Variable *atoms* specifies
        which atoms or rather which part of the molecule needs to be moved for
        this specific bond. The given length is going to be the new bond length,
        not by how much the bond length is changed.

        The move vector is determined automatically by the given length and atom bond.
        This vector can also be given manually with no regards to length,
        by setting the vector *vec*.

        Parameters
        ----------
        bond : list
            List of two atom ids of the bond to be adjusted
        atoms : int, list
            List of atoms that need to be moved by changing the bond length (can also be one id)
        length : float
            New bond length
        vec : list
            Set this vector, to manually move the atoms

        Examples
        --------
        >>> partMove([0,1],[1,2,3],0.5)
        """
        # Initialize
        data = self.__data
        dim  = self.__dim

        # Create temporary molecule
        if isinstance(atoms,int): atoms = [atoms]
        temp = self.__temp(atoms)

        # Set length
        length = abs(length - self.bond(bond))

        # Set vector
        if vec==None:
            vec = self.__vector(bond[0],bond[1])
        vec = [v*length for v in self.__unit(vec)]

        # Move molecule
        temp.translate(vec)

        # Update positions
        self.__update(temp.getData(),atoms)

    def partRotate(self,bond,atoms,angle,zero):
        """Rotate a set of specified atoms around a given bond as the rotation
        axis. First however the system needs to be set to zero. Therefore the
        atom id to define the new coordinate system has to be given for the set
        of specified atoms. Normally this is the atoms that lies at on an end
        of the given bond axis.

        Parameters
        ----------

        bond : list
            List of two atom ids of the bond to be set as an axis
        atoms : int, list
            List of atoms to be rotated (can also be one id)
        angle : float
            Rotation angle
        zero : int
            Atom id to define zero point of the new coordinate system

        Examples
        --------
        >>> partRotate([0,1],[1,2,3],90,0)
        """
        # Initialize
        data = self.__data
        dim  = self.__dim

        # Create temporary molecule
        self.move(zero,[0,0,0])
        if isinstance(atoms,int): atoms = [atoms]
        temp = self.__temp(atoms)

        # Rotate molecule
        temp.rotate([self.pos(bond[0]),self.pos(bond[1])],angle)

        # Update positions
        self.__update(temp.getData(),atoms)

    def partAngle(self,bondA,bondB,atoms,angle,zero):
        """Change the bond angle of two bond vectors. Variable *atoms* specifies
        which atoms or rather which part of the molecule needs rotated to change
        the specified bond angle. First however the system needs to be set to
        zero. Therefore the atom id to define the new coordinate system has to be
        given for the set of specified atoms. Normally this is the atom that
        touches the angle.

        The rotation axis is determined by creating the cross product
        of the two bond vectors. Thus getting the normal vector of a surface
        that contains both bond vectors.

        Parameters
        ----------
        bondA : list
            First bond vector given as a list of two atom ids
        bondB : list
            Second bond vector given as a list of two atom ids
        atoms : int, list
            List of atoms to be rotated (can also be one id)
        angle : float
            Rotation angle
        zero : int
            Atom id to define zero point of the new coordinate system

        Examples
        --------
        >>> partAngle([0,1],[1,2],[1,2,3],90,1)
        """
        # Initialize
        data = self.__data
        dim  = self.__dim

        # Create temporary molecule
        self.move(zero,[0,0,0])
        if isinstance(atoms,int): atoms = [atoms]
        temp = self.__temp(atoms)

        # Rotate molecule around normal vector
        if   len(bondA)==2   and len(bondB)==2:
            vec = self.__cross(self.__vector(bondA[0],bondA[1]),self.__vector(bondB[0],bondB[1]))
        elif len(bondA)==dim and len(bondB)==dim:
            vec = self.__cross(bondA,bondB)
        else:
            print("Wrong bond input...")
            return
        temp.rotate(vec,angle)

        # Update positions
        self.__update(temp.getData(),atoms)

    def lengthAngle(self,bond,bondA,bondB,atoms,zero,length,angle=[0,45],grid=0.001,isSilent=True,isNegative=False):
        """If the bond length of an exemplary circular molecule is correlated
        to an angle, then the bond length has to be changed in dependence to
        the bond angle.

        Variable *atoms* specifies which atoms or rather which part of the
        molecule needs moved and rotated to change the specified bond length.

        First however the system needs to be set to zero. Therefore the atom id
        to define the new coordinate system has to be given for the set of
        specified atoms. Normally this is the atom that touches the angle.

        The rotation axis is determined by creating the cross product
        of the two bond vectors. Thus getting the normal vector of a surface
        that contains both bond vectors.

        The new length is determined, by rotating the specified molecule part in
        steps and calculating the resulting length error. If the error is zero,
        then the rotation steps stop.

        Parameters
        ----------
        bond : list
            Bond of which the length should be changed
        bondA : list
            First bond vector given as a list of two atom ids
        bondB : list
            Second bond vector given as a list of two atom ids
        atoms : int, list
            List of atoms to be rotated (can also be one id)
        zero : int
            Atom id to define zero point of the new coordinate system
        length : float
            New bond length
        angle : list
            Angle range to be tested
        grid : float
            angle steps fineness
        isSilent : bool
            True to suppress error-value messages
        isNegative : bool
            True for a negative rotation

        Examples
        --------
        >>> lengthAngle([0,19],[8,9],[9,10],[10,11,12],self.pos(9))

        """
        angles = [x*grid for x in range(int(angle[0]/grid),int(angle[1]/grid)+1)]
        for a in angles:
            a = -a if isNegative else a
            self.partAngle(bondA,bondB,atoms,a,zero)

            error = round(abs(self.bond(bond)-length),3)

            if not isSilent:
                print("Error = "+str(error))

            if error==0:
                break


    #############################
    # Public Methods - Calulate #
    #############################
    def angle(self,vec1,vec2,isDeg=True):
        """Calculate the angle using function :func:`_angle`.

        Parameters
        ----------
        vec1 : list
            First vector
        vec2 : list
            Second vector
        isDeg : bool
            True if the output should be in degree

        Returns
        -------
        angle : float
            Angle value
        """
        return self.__angle(vec1,vec2,isDeg)

    def bond(self,atoms):
        """Calculate the bond length of a given bond (list of two atom ids).

        Parameters
        ----------
        atoms : list
            Bond atoms

        Returns
        -------
        bond : float
            Bond length
        """
        return self.__length(self.__vector(atoms[0],atoms[1]))

    def pos(self,atom):
        """Get the position of an atom as a vector.

        Parameters
        ----------
        atom : int
            Atom id

        Returns
        -------
        pos : list
            Position vector of the specified atom
        """
        pos = []
        for i in range(self.__dim):
            pos.append(self.__data[i][atom])

        return pos


    #########################
    # Public Methods - Edit #
    #########################
    def add(self,name,pos,bond=None,r=0,theta=0,phi=0,isDeg=True,bondType="s"):
        """Add a new atom in polar coordinates. The *pos* input is either
        an atom id, that determines is the bond-start,
        or a vector for a specific position.

        Bond has to be given, if the polar coordinates are dependet on the bond
        vector as the basic axis. The coordinate system is then transformed
        to the new bond axis. If set to None, then the given coordinates are
        thought to be dependent on the basic axes.

        If the given position is an atom id, then a bond is automatically added
        to the bond matrix.

        Parameters
        ----------
        name : str
            Atom type
        pos : int, list
            Position of the atom
        bond : None, list
            Bond axis
        r : float
            Bond length
        theta : float
            Azimuthal angle
        phi : float
            Polar angle
        isDeg : bool
            True if the input of the angles in degree
        bondType : str
            Bond type "s"-single, "d"-double

        Examples
        --------
        >>> add("C",[0,0,0])
        >>> add("C",0,r=0.153,theta=-135)
        >>> add("C",1,[0,1],r=0.153,theta= 135)
        """
        # Initialize
        data = self.__data

        # Angles
        phi   *= math.pi/180 if isDeg else 1
        theta *= math.pi/180 if isDeg else 1

        # Transform sphercial to cartesian coordinates
        x     = r*math.sin(theta)*math.cos(phi)
        y     = r*math.sin(theta)*math.sin(phi)
        z     = r*math.cos(theta)
        coord = [x,y,z]

        # Bond vector
        if bond==None:
            vec = self.__axis("z")
        else:
            vec = self.__vector(bond[0],bond[1])

        # Calculate angles for rotation
        phi   = self.__anglePolar(vec)
        theta = self.__angleAzi  (vec)

        # Rotate towards axis
        tRot = self.__cross([0,0,1],vec)
        if sum(tRot)==0: tRot  = "y"

        coord = self.__rotate(coord,tRot,theta,isDeg=False)
        #coord = self.__rotate(coord,"z",phi,  isDeg=False)
        #coord = [round(i,4) for i in coord]

        # Process position input
        b = None
        if isinstance(pos,int):
            b   = pos
            pos = self.pos(pos)


        # Add coordinates
        for i in range(self.__dim):
            data[i].append(pos[i]+coord[i])

        # Add name
        data[self.__dim].append(name)

    # Delete an atom
    def delete(self,atom):
        """Delete specified atom from the molecule.
        The input can also be a list of atom ids.

        Parameters
        ----------
        atom :  int, list
            Atom id or list to be deleted
        """
        # Initialize
        data = self.__data

        # Process input
        if isinstance(atom,str):
            if atom == "last":
                atom = [len(data[0])-1]
        elif not isinstance(atom,list):
            atom = [atom]
        atom = sorted(atom,reverse=True)

        # Remove line
        for dat in data:
            for a in atom:
                dat.pop(a)

    def overlap(self,error=0.005,isPrint=False):
        """Check if atoms are overlapping (as in duplicates) and delete these.

        Parameters
        ----------
        error : float
            Error of the overlap search
        isPrint : bool
            True to print the list of overlapping atoms
        """
        # Initialize
        data = self.__data
        dim  = self.__dim

        # Find duplicate atoms in reverse
        for i in list(reversed(range(self.getNum()))):
            # Get positions
            atomA = self.pos(i)
            for j in range(self.getNum()):
                if not j==i:
                    atomB  = self.pos(j)

                    # Check if identical
                    sumDel = 0
                    a = []
                    for k in range(dim):
                        if abs(atomA[k]-atomB[k]) < error:
                            a.append(abs(atomA[k]-atomB[k]))
                            sumDel += 1

                    # Delete atom
                    if sumDel == dim:
                        # print(i,j,a)
                        self.delete(i)
                        if isPrint:
                            print(i)
                        break

    def bondLength(self,inpA,inpB):
        """Return the bond length of a specified bond. The two inputs can either
        be atom indices or to vectorial positions.

        Parameters
        ----------
        inpA : int, list
            Either an atom id or a position vector
        inpB : int, list
            Either an atom id or a position vector

        Returns
        -------
        length : float
            Bond length

        Examples
        --------
        >>> bondLength(0,1)
        >>> bondLength([1,0,0],[0,0,0])
        """
        return self.__length(self.__vector(inpA,inpB))

    def bondVec(self,inpA,inpB):
        """Return the bond vector of a given bond. the Two inputs can either
        be atom indices or to vectorial positions.

        Parameters
        ----------
        inpA : int, list
            Either an atom id or a position vector
        inpB : int, list
            Either an atom id or a position vector

        Returns
        -------
        bond : list
            Bond vector
        """
        return self.__vector(inpA,inpB)

    def changeType(self,atom,atomType):
        """Change the atom type of a specified atom.

        Parameters
        ----------
        atom : int
            Atom id
        atomType : str
            New atomtype
        """
        self.__data[self.__dim][atom] = atomType

    def getMin(self,inp):
        """Return minimal coordinate value for the specified
        dimension (start at zero).

        Parameters
        ----------
        inp : inst
            Axis id

        Returns
        -------
        minimum : float
            Minimal coordinate
        """
        # Initialize
        data = self.__data

        # Check minimum
        minimum = 0
        for i in range(len(data[inp])):
            if minimum>data[inp][i]:
                minimum = data[inp][i]

        return minimum


    ##################
    # Setter Methods #
    ##################
    def setName(self,name):
        """Set the molecule name, short name and link based on that name.

        Parameters
        ----------
        name : str
            Molecule name
        """
        self.__name = name

    def setLink(self,link=None):
        """Set the molecule folder link. If input is None, set the link automatically.

        Parameters
        ----------
        link : str, None
            Molecule folder link
        """
        self.__link = config.load("system")["home"]+config.load("ident")["links"]["mols"]["out"]+self.getName()+"/" if link is None else link

    def setShort(self,short=None):
        """Set the molecule short name. If input is None, set the short name automatically.

        Parameters
        ----------
        short : str, None
            Molecule short name
        """
        self.__short = config.load("mols")[self.getName()] if short is None else short

    def setWrite(self,write):
        """Set the molecule list for writing the structure file

        Parameters
        ----------
        write : list
            List of molecule objects
        """
        self.__write = write

    def setBox(self,box):
        """Set the molecule box dimensions.

        Parameters
        ----------
        box : list
            Box dimension
        """
        self.__box = box

    def setCharge(self,charge):
        """Set the total charge of the molecule.

        Parameters
        ----------
        charge : float
            Total molecule charge
        """
        self.__charge = charge

    def setMasses(self,masses=None):
        """Set the molar masses of the atoms.

        Parameters
        ----------
        masses : list
            List of molar masses in :math:`\\frac g{mol}`
        """
        if masses is not None:
            self.__masses = masses
        else:
            self.__masses = []
            for atom in self.__data[self.__dim]:
                self.__masses.append(self.__db.getMass(atom))

    def setMass(self,mass=None):
        """Set the molar mass of the molecule.

        Parameters
        ----------
        mass : float
            Molar mass in :math:`\\frac g{mol}`
        """
        if mass is not None:
            self.__mass = mass
        else:
            self.__mass = 0
            for mass in self.getMasses():
                self.__mass += mass

    def setCOM(self,com=None):
        """Set the center of mass coordinates of the molecule.

        Parameters
        ----------
        com : list
            Center of mass coordinates
        """
        if com is not None:
            self.__com = com
        else:
            self.__com = []
            masses    = self.getMasses()
            for i in range(self.__dim):
                self.__com.append(0)
                for j in range(len(self.__data[i])):
                    self.__com[i] += masses[j]*self.__data[i][j]

            self.__com = [x/sum(masses) for x in self.__com]


    ##################
    # Getter Methods #
    ##################
    def getData(self):
        """Return the data matrix.

        Returns
        -------
        data : list
            Data matrix
        """
        return self.__data

    def getName(self):
        """Return the molecule name.

        Returns
        -------
        name : str
            Molecule name
        """
        return self.__name

    def getLink(self):
        """Return the molecule folder link.

        Returns
        -------
        link : str
            Molecule folder link
        """
        if self.__link is None: self.setLink()

        return self.__link

    def getShort(self):
        """Return the molecule short name.

        Returns
        -------
        short : str
            Molecule short name
        """
        if self.__short is None: self.setShort()

        return self.__short

    def getWrite(self):
        """Return the global molecule list for writing the structure file.

        Returns
        -------
        write : list
            List of all molecules
        """
        return self.__write

    def getBox(self):
        """Calculate and return the box size of the current molecule.

        Returns
        -------
        box : list
            Box dimensions
        """
        return self.__boxSize()

    def getBoxC(self):
        """Return the current box size.

        Returns
        -------
        box : list
            Box dimensions
        """
        return self.__box

    def getNum(self):
        """Return the number of atoms.

        Returns
        -------
        num : int
            Number of atoms
        """
        return len(self.__data[0])

    def getType(self,atom):
        """Return the atomtype of the given atom id.

        Parameters
        ----------
        atom : int
            Atom id

        Returns
        -------
        box : list
            Box dimensions
        """
        return self.__data[self.__dim][atom]

    def getCharge(self):
        """Return the total charge of the molecule.

        Returns
        -------
        charge : float
            Total charge
        """
        return self.__charge

    def getMasses(self):
        """Return a list of masses of the atoms.

        Returns
        -------
        masses : list
            Masses in :math:`\\frac g{mol}`
        """
        if self.__masses is None: self.setMasses()
        return self.__masses

    def getMass(self):
        """Return the molar mass of the molecule.

        Returns
        -------
        mass : float
            Molar mass in :math:`\\frac g{mol}`
        """
        if self.__mass is None: self.setMass()
        return self.__mass

    def getCOM(self):
        """Return the center of mass of the molecule.

        Returns
        -------
        com : list
            center of mass coordinates
        """
        if self.__com is None: self.setCOM()
        return self.__com
