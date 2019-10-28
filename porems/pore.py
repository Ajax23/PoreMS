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
from porems.verlet   import Verlet
from porems.bonding  import Bonding


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

    For an **exemplary** creation consider class extension :class:`pores.mols.epoxiPore.EpoxiPore`.

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
    def __init__(self,size,diam,drill,res=1,vs=0.4,is_pbc=True,is_time=False):
        # Call super class
        super(Pore, self).__init__()

        # Initialize
        self._dim     = config.load("prefs")["dim"]
        self._size    = size if isinstance(size,list) else [size for s in range(self._dim)]
        self._diam    = diam-{"x":0.55,"y":0.5,"z":0.5}[drill]
        self._is_pbc  = is_pbc
        self._is_time = is_time
        self._res     = res
        self._vs      = vs
        self._drill   = drill

        # Define bond lengths
        self._siGrid  = self._db.getBond("Si","O")
        self._siLol   = self._db.getBond("Si","OH")
        self._siLox   = self._db.getBond("Si","OH")
        self._oh      = self._db.getBond("O", "H")

        # Set proximity range
        self._oGrid   = 0.300

        # Define silica block information
        self._repeat  = [0.506,0.877,1.240]
        self._gap     = [0.126,0.073,0.155]
        self._numRep  = [round(self._size[i]/self._repeat[i])         for i in range(self._dim)]
        self._size    = [self._repeat[i]*self._numRep[i]+self._gap[i] for i in range(self._dim)]

        # Define charges
        self._qSi     =  1.28
        self._qO      = -0.64
        self._qOH     = -0.74+0.42

        # Define lists
        self._site    = []
        self._grid    = []

        # Time management
        self._tTot    = 0

        # Create silica pore
        t = m.tic()
        self._block(0,self._foundation())              # Build block
        self._orientation()                            # Rotate drill axis
        self.translate(self._gap)                      # Translate gap
        self._center  = self._focal()                  # Find focal point
        self._tTot   += m.toc(t,"Build  ",is_time)

        # Verlet and bonding
        t = m.tic()
        self._verlet   = Verlet(self,vs,is_pbc)         # Create verlet boxes
        self._bonding  = Bonding(self._verlet)         # Create bond matrix
        self._tTot    += m.toc(t,"Matrix ",is_time)

        # Prepare pore
        t = m.tic()
        self._box   = self.getBox()                    # Set box
        self._bonding.attach()                         # Preperare sides
        self._bonding.drill(self._center,self._diam)   # Drill pore
        self._bonding.prepare()                        # Prepare pore surface
        self._tTot += m.toc(t,"Prepare",is_time)

        # Find binding sites
        t = m.tic()
        self.zero()
        self._box = self.getBox()                      # Reset box size
        self._bind()                                   # Create site array
        self._proxi()                                  # Find sites in proximity
        self._diam = self._diameter()                  # Calculate new diameter
        self._tTot += m.toc(t,"Binding",is_time)


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
        dim   = self._dim
        dist  = self._siGrid
        angle = 2*math.atan(math.sqrt(2))*180/math.pi

        # Create pore
        h00 = Molecule()
        h00.add("Si",[0,0,0])
        h00.add("O" , 0,r= dist,theta=angle,phi= 60)
        h00.add("Si", 1,bond=[0,1],r=dist)
        h00.add("O" , 2,r= dist,theta=180,  phi=  0)
        h00.add("Si", 3,bond=[2,3],r=dist)
        h00.add("O" , 4,r= dist,theta=angle,phi=300)
        h00.add("Si", 5,bond=[4,5],r=dist)
        h00.add("O" , 6,r=-dist,theta=angle,phi= 60)
        h00.add("Si", 7,bond=[6,7],r=dist)
        h00.add("O" , 8,r=-dist,theta=180,  phi=  0)
        h00.add("Si", 9,bond=[8,9],r=dist)
        h00.add("O" ,10,r=-dist,theta=angle,phi=300)

        h00.rotate("y", 90)
        h00.rotate("z", 90)
        h00.rotate("y",180)
        h00.rotate("x",h00.angle(h00.bondVec(8,6),self._axis("y")))
        h00.rotate("x",-90)

        h01 = copy.deepcopy(h00)
        h01.add("O", 2,r=dist)
        h01.add("O", 6,r=dist)
        h01.add("O",10,r=dist)

        h02 = copy.deepcopy(h00)
        h02.move(0,h01.pos(12))
        h02.translate([0,0,dist])
        h02.delete([1,2,3,4,5,6,7])

        h03 = copy.deepcopy(h00)
        h03.move(0,h01.pos(14))
        h03.translate([0,0,dist])
        h03.delete([2,3,4,5,6,7,8,9,10,11])

        hexa = Molecule([h01,h02,h03])

        # Calculate xrepeat
        z  = [hexa]
        for i in range(10): z.append(copy.deepcopy(hexa))

        z[1].rotate("x",180)
        z[1].move(18,z[0].pos(18))
        z[1].translate([0,0,dist*2])
        z[1].add("O",z[0].pos(18),r=dist)

        z[2].rotate("x",180)
        z[2].move(0,z[0].pos(18))

        z[3].move(4,z[0].pos(15))

        z[4].rotate("x",180)
        z[4].move(16,z[0].pos(18))
        z[4].delete([21,19,15,20,14,12,10,2,11,1,0])

        z[5].move(6,z[1].pos(16))
        z[5].delete([21,19,15,20,14,12,10,2,11,1,0])

        # Create repeat block
        block = Molecule(z)
        block.overlap()

        # Determine number of needed blocks
        block.zero()

        # Delete overlapping atoms
        block.delete([4,3,2,12,15,45,46,55,57,63,26,25,24,34,37,61]) # x-Axis
        block.delete([50,48,36,40,41])                               # y-Axis
        block.delete([45,20,19,21,22,23,24,25,17,18])                # z-Axis

        # Move to zero
        block.zero()

        return block

    def _block(self,dim,block):
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

            for i in range(self._numRep[dim]):
                temp = copy.deepcopy(block)
                vec  = [(i+1)*self._repeat[dim] if j==dim else 0 for j in range(self._dim)]

                temp.translate(vec)

                p.append(temp)

            self._block(dim+1,Molecule(p))
        else:
            self._append(block)

    def _orientation(self):
        """Rotate pore orientation, so that the specified drill axis becomes
        the z-axis.
        """
        # Initialize
        drill = self._drill
        gap   = self._gap
        size  = self._size

        # Rotate pore
        if   drill=="x": self.rotate("y",90)
        elif drill=="y": self.rotate("x",90)

        # Set zero
        self.zero()

        # Update gap and size lists
        if   drill=="x":
            self._gap  = [gap [2],gap [1],gap [0]]
            self._size = [size[2],size[1],size[0]]
        elif drill=="y":
            self._gap  = [gap [0],gap [2],gap [1]]
            self._size = [size[0],size[2],size[1]]


    ###################################
    # Private Methods - Binding sites #
    ###################################
    def _bind(self):
        """Create binding site matrix. This matrix :math:`\\boldsymbol{B}`
        is mostly generated by the bonding class :class:`pores.build.bonding.Bonding`

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
        osi  = self._bonding.site()

        # Create bonding site list
        for os in osi:
            site.append([os[0],os[1],[],0,os[2],os[3],-1])

        self._site.extend(site)

    def _proxi(self):
        """Using verlet lists, search for binding sites that are in proximity,
        by looking for oxygene atoms that are near each other.

        Fill the result in the binding site matrix in the proximity entry.
        """
        # Initialize
        site    = self._site

        # Create local verlet list
        atoms  = [o[0] for o in site]
        verlet = Verlet(self._temp(atoms),self._vs,self._is_pbc)

        # Find oxygenes in proximity
        boxList = [i for i in range(len(verlet.getBox()[1]))]
        oo      = verlet.findParallel(boxList,["O","O"],self._oGrid,10e-2)
        ooCol   = m.column(oo)

        # Append proxi list to sites
        for i in range(len(site)):
            os = site[i]
            o  = oo [ooCol[0].index(i)][1] if i in ooCol[0] else []

            os[2].extend(o)

    def _random(self,siteType,rate,inp,counter,isProxi=False):
        """Get a random binding site list for a given percentage of the available binding sites.

        If a site is available (state is 0), it is added to the output-list and
        the state entry for this site is set to 1.
        Then the states of all binding sites in the proximity list are set to 2.

        In case the random site is not available, this run is repeated. A break-counter
        prevents the occurance of an endless loop.

        If *isProxi* is set to True, then binding sites nearest to each other
        are determined. the return value then consists of lists containing the
        ids of both sites. Geminal pairs are not allowed, since creating the
        topology is too complex. This however might be added in the future.

        Parameters
        ----------
        siteType : int
            Type of the sites inside-0, outside-1
        rate : float
            Rate of how many binding sites
        inp : str
            Input type of the rate
        counter : int
            Number of attempts before breaking the randoming loop
        isProxi : bool
            True to search for site pairs in proximity

        Returns
        -------
        rand : list
            List of pointers to random binding sites.
        """
        # Initialize
        site    = self._site
        temp    = m.column(site)
        siteLen = temp[4].count(siteType)
        idList  = [i for i in range(len(site))]

        # Calculate random element number
        if   inp=="num":
            ranNum  = rate
        elif inp=="percent":
            ranNum  = siteLen*rate
            ranNum /= 100
            ranNum  = math.floor(ranNum)

        # Get list of all partners
        if isProxi:
            siteMin = []
            for i in range(len(site)):
                length = [self._length(self._vector(site[i][0],site[j][0])) if not i==j and site[i][4]==site[j][4] else 100000 for j in range(len(site))]

                siteMin.append([length.index(min(length)),min(length)])

        # Get random elements from site list
        ranList = []
        count   = 0
        while len(ranList)<ranNum and count<counter:
            # Get random site index
            rand = random.choice(idList)

            # Check binding site
            isApp = (isProxi             and site[rand][4]==siteType      and
                     site[rand][3]==0    and site[siteMin[rand][0]][3]==0 and
                     site[rand][5]==None and site[siteMin[rand][0]][5]==None)

            isApp = site[rand][4]==siteType and site[rand][3]==0 if not isProxi else isApp

            # Append to return list
            if isApp:
                # Reset counter
                count = 0

                # Append
                ranList.append([rand,siteMin[rand][0]] if isProxi else rand)

                # Set site and proximity to used
                self._close([rand,siteMin[rand][0]] if isProxi else rand)
            else:
                # Add counter
                count += 1

        return ranList

    def _close(self,sites):
        """This function cloese the specified site.

        Parameters
        ----------
        sites : int, list
            identifiers of sites to be closed
        """
        # Initialize
        site = self._site

        # Process input
        sites = [sites] if isinstance(sites,int) else sites

        # Close sites
        for x in sites:
            site[x][3] = 1
            for prox in site[x][2]:
                if not site[prox][3]==1:
                    site[prox][3] = 2


    ##################################
    # Private Methods - Add Molecule #
    ##################################
    def _couple(self,mol,sio,siteType,rate,inp="percent",sites=None,counter=1000):
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
        siteType : int
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
        dim   = self._dim
        site  = self._site

        # Process user input
        sites = [sites] if isinstance(sites,int) else sites

        # Get random list or given list
        ranList = self._random(siteType,rate,inp,counter) if sites is None else sites

        # Molecule and center vector
        vecM   = mol.bondVec(sio[1],sio[2])
        center = [self._center[i] for i in range(dim-1)]

        # Add molcule
        writeL = [[],[]]
        for i in sorted(ranList,reverse=True):
            # Initialize
            o = site[i]

            # Copy molecule
            temp = copy.deepcopy(mol)

            # Calculate osi vector
            vecO = self._vector(o[1],o[0])
            vecC = center[:]
            vecC.append(self.pos(o[0])[dim-1])
            vecC = self._vector(self.pos(o[1]),vecC)

            # Random axis rotation
            temp.rotate(vecM,random.randint(1,180))

            # Adjust molecule vector to osi
            temp.rotate(self._cross(vecM,vecO),self._angle(vecM,vecO))

            # Adjust molecule towards center
            if siteType==0:
                temp.rotate(self._cross(vecC,vecO),-self._angle(vecC,vecO))

            # Move molecule to position
            temp.move(sio[1],self.pos(o[0]))

            # Move silica
            temp.put(sio[0],self.pos(o[1]))

            # Check if geminal
            isGem = not o[5] == None

            name = temp.getName()+"g" if isGem else temp.getName()
            temp.setName(name)

            # Remove atoms
            self._bonding.remove([o[0],o[1]])

            # Add molecule to pore
            writeID = 1 if isGem else 0
            writeL[writeID].append([temp,i])

        # Add to write
        for write in writeL:
            for i in range(len(write)):
                site[write[i][1]][6] = len(self._write)
                self._write.append(write[i][0])

    def _coupleProxi(self,mol,sio,orient,rate,inp="percent",counter=1000):
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
        dim  = self._dim
        site = self._site

        # Get random list
        ranList = self._random(0,rate,inp,counter,isProxi=True)

        # Molecule and center vector
        center = [self._center[i] for i in range(dim-1)]

        # Add molcule
        for i in sorted(ranList,reverse=True):
            # Initialize
            s = [site[i[0]],site[i[1]]]

            # Copy molecule
            temp = copy.deepcopy(mol)

            # Rotate towards binding site
            vecMO = temp.bondVec(sio[0][1],sio[1][1]) # Molecule oxygenes
            vecSO = self._vector(s[0][0],s[1][0])     # Binding site oxygenes
            vecM  = temp.bondVec(orient[0],orient[1]) # Molecule orientation

            temp.rotate(vecM,self._angle(vecMO,vecSO))

            # Rotate towards central axis
            centO = [(self.pos(s[0][0])[x]+self.pos(s[1][0])[x])/2 for x in range(dim)]
            centP = center[:]+[centO[-1]]

            vecM  = temp.bondVec(orient[0],orient[1]) # Molecule orientation
            vecC  = self._vector(centO,centP)         # Binding site center

            temp.rotate(self._cross(vecM,vecC),self._angle(vecM,vecC))

            # Put oxygenes ontop of each other
            temp.move(sio[0][1],self.pos(s[0][0]))

            # Move to center
            centM = [(temp.pos(sio[0][1])[x]+temp.pos(sio[1][1])[x])/2 for x in range(dim)]
            temp.translate(self._vector(centM,centO))

            # Move silica atoms
            temp.put(sio[0][0],self.pos(s[0][1]))
            temp.put(sio[1][0],self.pos(s[1][1]))

            # Remove atoms
            self._bonding.remove([s[0][0],s[0][1],s[1][0],s[1][1]])

            # Add to write
            self._write.append(temp)

    def _special(self,mol,sio,num,symmetry):
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
        diam   = self._diam
        site   = self._site
        center = self._center

        # Process input
        if symmetry not in ["random","point","mirror"]:
            print("Symmetry type not supported...")
            return

        elif not symmetry=="random" and not num==0:
            # Calculate placment points
            posL   = []
            length = center[2]*2
            dist   = length/num
            start  = dist/2

            for i in range(num):
                if   symmetry=="point":
                    coeff = -1 if i%2==0 else 1
                elif symmetry=="mirror":
                    coeff = 1

                x = center[0]+coeff*diam/2
                y = center[1]
                z = start+dist*i

                posL.append([x,y,z])

            # Find nearest site
            siteMin = [[0,None] for i in range(num)]
            for i in range(len(site)):
                length = [self._length(self._vector(self.pos(site[i][0]),pos)) for pos in posL]

                for j in range(num):
                    if siteMin[j][1] is None or length[j]<siteMin[j][1]:
                        siteMin[j] = [i,length[j]]

            sites = [x[0] for x in siteMin]

            self._close(sites)

        else:
            sites = None

        self._couple(mol,sio,0,num,inp="num",sites=sites)


    ################################
    # Private Methods - Final Mols #
    ################################
    def _silanol(self,siteRange=None):
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
        writeL : list
            List of Molecule objects and ids.
        """
        # Initialize
        site   = self._site
        dim    = self._dim
        box    = self._box
        center = [self._center[i] for i in range(dim-1)]
        writeL = [[],[]]

        # Set siterange
        siteR = [0,len(site)] if siteRange is None else siteRange

        # Define temporary molecule object
        tempMol = Molecule()

        # Change bond length and add hydrogen
        for i in range(siteR[0],siteR[1]):
            if not site[i][3]==1:
                # Initialize
                o    = site[i]
                o[3] = 1
                temp = copy.deepcopy(tempMol)

                # Change SiO bond length
                if self.bondLength(o[0],o[1])<box[2]/2:
                    self.partMove([o[1],o[0]],o[0],self._siLol)

                # Inside pore
                if o[4]==0:
                    temp.add("O",[0,0,0])
                    temp.add("H",0,r=self._oh)

                    posC = center[:]
                    posC.append(self.pos(o[0])[dim-1])

                    vecM = temp.bondVec(0,1)
                    vecC = self._vector(self.pos(o[1]),posC)
                    vecO = self._vector(o[1],o[0])

                    temp.rotate(self._cross(vecM,vecO), self._angle(vecM,vecO))
                    temp.rotate(self._cross(vecC,vecO),-self._angle(vecC,vecO))

                    temp.move(0,self.pos(o[0]))
                    temp.add("Si",self.pos(o[1]))

                    # Sort
                    tempPosO = temp.pos(0)
                    tempPosH = temp.pos(1)
                    temp.delete([0,1])
                    temp.add("O",tempPosO)
                    temp.add("H",tempPosH)

                # Outside pore
                else:
                    angle = -180 if self.pos(o[0])[2]<box[2]/2 else 0

                    temp.add("Si",self.pos(o[1]))
                    temp.add("O", self.pos(o[0]))
                    temp.add("H",1,r=self._oh,theta=angle)

                # Check if geminal
                isGem = not o[5]==None

                # Add Charge
                temp.setCharge(self._qSi+self._qOH)

                # Add to write
                sl = "slg" if isGem else "sl"
                temp.setName(sl)

                writeID = 1 if isGem else 0
                writeL[writeID].append([temp,i])

        # Add to write and add site - molecule connection
        if siteRange is None:
            for write in writeL:
                for w in write:
                    self._bonding.remove([site[w[1]][1],site[w[1]][0]])
                    site[w[1]][6] = len(self._write)
                    self._write.append(w[0])
        else:
            return writeL

    def _silanolParallel(self):
        """Parallelized function :func:`_silanol`.
        """
        # Initialize
        np    = config.load("system")["np"]
        sites = self._site

        # Split site list
        siteLen  = math.floor(len(sites)/np)
        siteList = []
        for i in range(np):
            if i==np-1:
                siteList.append([siteLen*i,len(sites)])
            else:
                siteList.append([siteLen*i,siteLen*(i+1)])

        # Paralellize
        pool    = mp.Pool(processes=np)
        results = pool.map_async(self._silanol,siteList)
        pool.close()

        writeL = [[],[]]
        for result in results.get():
            writeL[0].extend(result[0])
            writeL[1].extend(result[1])

        # Add to write and add site - molecule connection
        for write in writeL:
            for w in write:
                self._bonding.remove([sites[w[1]][1],sites[w[1]][0]])
                sites[w[1]][6] = len(self._write)
                self._write.append(w[0])

    # Transform two SiO to siloxan bridges
    def _siloxan(self,rate,inp="percent",counter=1000):
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
        dim  = self._dim

        # Exit if rate is zero
        if rate==0: return

        # Get random list
        ranList = self._random(0,rate,inp,counter,isProxi=True)

        # Molecule and center vector
        center = [self._center[i] for i in range(dim-1)]

        # Define temporary molecule object
        tempMol = Molecule()

        # Add molcule
        for i in sorted(ranList,reverse=True):
            # Initialize
            s = [site[i[0]],site[i[1]]]

            # Get center position
            centO = [(self.pos(s[0][0])[x]+self.pos(s[1][0])[x])/2 for x in range(dim)]
            centS = [(self.pos(s[0][1])[x]+self.pos(s[1][1])[x])/2 for x in range(dim)]

            # Create molecule
            temp = copy.deepcopy(tempMol)
            temp.setName("slx")
            temp.add("OM",centO)
            temp.translate(self._vector(centS,centO))
            temp.add("O",self.pos(s[0][0]))
            temp.add("O",self.pos(s[1][0]))

            # Rotate towards central axis
            centP = center[:]+[centO[-1]]

            vecM  = self._vector(centO,temp.pos(0))
            vecC  = self._vector(centO,centP)

            temp.rotate(self._cross(vecM,vecC),self._angle(vecM,vecC))

            # Put oxygenes ontop of each other
            temp.move(1,self.pos(s[0][0]))

            # Remove temporary oxygen atoms
            temp.delete([1,2])

            # Remove atoms
            self._bonding.remove([s[0][0],s[1][0]])

            # Add to write
            self._write.append(temp)


    #################################
    # Private Methods - Final Edits #
    #################################
    def _connect(self):
        """Concatenate two geminal molecules to one molecule.
        """
        # Initialize
        site  = self._site
        write = self._write
        dim   = self._dim
        gemi  = []

        # Check for geminal sites
        for i in range(len(site)):
            if site[i][5] is not None:
                gem  = [site[i][6],site[site[i][5]][6]]
                gemR = [site[site[i][5]][6],site[i][6]]

                if not gemR in gemi:
                    gemi.append(gem)

        # Append atoms
        pop = []
        for g in gemi:
            idA  = 0 if write[g[1]].getName()=="slg" else 1
            idB  = abs(idA-1)

            molA = write[g[idA]]
            molB = write[g[idB]]

            num  = molB.getNum()

            for i in range(num):
                if not molB.getType(i)=="Si":
                    molA.add(molB.getType(i),molB.pos(i))

            pop.append(g[idB])

            # Add charge
            write[g[idA]].setCharge(write[g[idA]].getCharge()+self._qOH)

        # Delete from write list
        for d in sorted(pop,reverse=True):
            write.pop(d)

    def _objectify(self):
        """Move all remaining grid silica and oxygen atoms to individual molecules
        for writing the structure file and add these molecules to the global
        molecule list after sorting.
        """
        # Initialize
        data   = self._data
        dim    = self._dim
        writeL = [[],[]]

        # Define temporary molecule object
        tempMol = Molecule()

        # Create OM and SI mols
        for i in range(len(data[0])):
            temp  = copy.deepcopy(tempMol)

            # Check if oxygene
            isOxy = self.getType(i)=="O"

            # Set temporary vars
            atom    = "OM"     if isOxy else "SI"
            name    = "om"     if isOxy else "si"
            charge  = self._qO if isOxy else self._qSi
            writeID = 0        if isOxy else 1

            # Testing
            if self.getType(i)=="R":
                atom = "R"
                name = "test"

            # Creat unique object and add to write
            temp.add(atom,self.pos(i))
            temp.setName(name)
            temp.setCharge(charge)
            writeL[writeID].append(temp)

        # Add mols to write
        self._write.pop(0)

        for write in writeL:
            self._write = write + self._write

    def _sort(self):
        """Sort molecules in order to prevent multiple molecule appearances in
        the structure file.
        """
        # Initialize
        write    = self._write
        writeNew = []

        # Get unique mols
        uniqueMols = []
        for mol in write:
            if mol.getName() not in uniqueMols:
                uniqueMols.append(mol.getName())

        # Sort molecules
        for molName in uniqueMols:
            for mol in write:
                if mol.getName()==molName:
                    writeNew.append(mol)

        self._write = writeNew

    def _position(self):
        """Tranlate the pore into position, and consider the periodic repeat gap.
        Hereby the coordinates are moved by halve the gap distance and box box is
        extended by the other halve. Additionally, to prevent molecule overlap and
        perodic movement, the drill side is extended and translated by a fixed value,
        thus creating reservoirs.
        """
        # Initialize
        dim   = self._dim
        box   = self._box
        gap   = self._gap
        write = self._write

        # Get vector to coordinate zero
        vec = Molecule([x for x in write if x.getName() in ["si","om"]]).zero()

        # Position pore
        for x in write:
            x.translate([vec[i]+gap[i]/2 for i in range(dim)])
            x.translate([self._res if i==2 else 0 for i in range(dim)])

        # Extend box
        coord = Molecule([x for x in write if x.getName() in ["si","om","ox"]]).getBox()

        for i in range(dim): box[i] = coord[i] + gap[i]/2
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
        verlet = Verlet(temp,self._vs,True)

        # Check silica and oxygene
        si = verlet.findParallel(None,["Si","Si"],0,10e-3)
        oo = verlet.findParallel(None,["O" ,"O"], 0,10e-3)

        # Output
        for s in si:
            if len(s[1])>0: print(s)

        for o in oo:
            if len(o[1])>0: print(o)

    def _excess(self,isRand=True):
        """Adjust grid atoms charges to compensate for excess charges on the pore.

        Parameters
        ----------
        isRand : bool
            True to randomly distribute rounding error charge on block oxygen atoms
        """
        # Initialize
        mols  = self._write
        grid  = [x.getName() for x in self._grid]+["sl"]
        grid += [x+"g" for x in grid]

        # Get excess charge
        charge = 0.0
        numO   = numSi = numG = 0
        roundO = 6
        for mol in mols:
            charge += mol.getCharge()

            if   mol.getName()=="om":   numO  += 1
            elif mol.getName()=="si":   numSi += 1
            elif mol.getName() in grid: numG  += 1

        # Calculate adjustment
        avg = charge/(numO/2+numSi+numG)
        qO  = round(self._qO -avg/2,roundO)
        qSi = round(self._qSi-avg,  roundO)

        # Distribute charge
        for mol in mols:
            if   mol.getName()=="om":   mol.setCharge(qO)
            elif mol.getName()=="si":   mol.setCharge(qSi)
            elif mol.getName() in grid: mol.setCharge(round(mol.getCharge()-self._qSi+qSi,roundO))

        # Calculate rounding error
        chargeErr = sum([mol.getCharge() for mol in mols])
        numOxy    = 100 if numO>100 else numO
        avgErr    = chargeErr/numOxy
        qOx       = round(qO-avgErr,roundO)

        # Randomly distribute balancing oxygens
        if isRand:
            # Get list of oxygenes
            oxy = {}
            for i in range(len(mols)):
                if mols[i].getName()=="om":
                    oxy[i] = mols[i]

            # Get random oxygenes and set name and charge
            oxyR = random.sample([x for x in oxy],numOxy)
            for o in oxyR:
                oxy[o].setName("ox")
                oxy[o].setCharge(qOx)

            # Push oxy after om
            minO    = min([x for x in oxy])
            maxO    = max([x for x in oxy])
            write   = []
            writeOx = []
            for i in range(minO): write.append(mols[i])
            for i in range(minO,maxO+1):
                if i not in oxyR: write.append(mols[i])
                else: writeOx.append(mols[i])
            write += writeOx
            for i in range(maxO+1,len(mols)): write.append(mols[i])

            self._write = write

        # Make global
        self._q = {"om":qO, "ox":qOx, "si":qSi}

        self.setCharge(sum([mol.getCharge() for mol in self._write]))


    ############################
    # Public Methods - Analyze #
    ############################
    def _allocation(self,site=None,isMol=True):
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
        isMol : bool
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
        site = self._site if site==None else site
        unit = 10/6.022 if isMol else 1

        # Calculate surfaces
        surface = {"pore":  size[2]*2*math.pi*(diam/2),
                   "drill": 2*(size[0]*size[1]-math.pi*(diam/2)**2)}

        # Get total number of sites
        count = {}
        count["pore"]  = {"tot" : sum([1 for x in site if x[4]==0]),
                          "mod":  sum([1 for x in site if x[4]==0 and x[3]==1])}
        count["drill"] = {"tot" : sum([1 for x in site if x[4]==1]),
                          "mod":  sum([1 for x in site if x[4]==1 and x[3]==1])}

        # Total allocation
        alloc = {}
        alloc["pore"]  = {"tot" : count["pore"]["tot"]/surface["pore"]*unit,
                          "mod":  count["pore"]["mod"]/surface["pore"]*unit,
                          "oh"  : (count["pore"]["tot"]-count["pore"]["mod"])/surface["pore"]*unit,
                          "rel" : count["pore"]["mod"]/count["pore"]["tot"]}
        alloc["drill"] = {"tot" : count["drill"]["tot"]/surface["drill"]*unit,
                          "mod":  count["drill"]["mod"]/surface["drill"]*unit,
                          "oh"  : (count["drill"]["tot"]-count["drill"]["mod"])/surface["drill"]*unit,
                          "rel" : count["drill"]["mod"]/count["drill"]["tot"]}

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
        site   = self._site
        center = self._center

        # Run through silica atoms and calculate distances
        rList = []
        for x in site:
            if x[4]==0:
                surf = x[1]
                cent = [center[0],center[1],self._data[2][surf]]
                rList.append(self._length(self._vector(self.pos(surf),cent)))

        # Calculate mean
        rBar = sum(rList)/len(rList)

        # Calculate roughness
        return math.sqrt(sum([(ri-rBar)**2 for ri in rList])/len(rList))

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
        site   = self._site
        center = self._center

        # Run through silica atoms and calculate distances
        rList = []
        for x in site:
            if x[4]==0:
                surf = x[1]
                cent = [center[0],center[1],self._data[2][surf]]
                rList.append(self._length(self._vector(self.pos(surf),cent)))

        # Calculate mean
        rBar = sum(rList)/len(rList)

        # Calculate roughness
        return 2*rBar


    #########################
    # Public Methods - Edit #
    #########################
    def special(self,mol,sio,num=2,symmetry="point"):
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
        self._special(mol,sio,num,symmetry)

    def couple(self,mol,sio,rate,intent="inside",inp="percent"):
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
        isList = isinstance(mol,list) and isinstance(sio,list) and isinstance(rate,list)
        if   intent=="inside":  intent = 0
        elif intent=="outside": intent = 1
        else:
            print("Wrong intent in coupling function...")
            return

        # Run coupling
        t = m.tic()

        if isList:
            if mol[0] is not None: self._couple(mol[0],sio[0],0,rate[0],inp=inp)
            if mol[1] is not None: self._couple(mol[1],sio[1],1,rate[1],inp=inp)
        else:
            self._couple(mol,sio,intent,rate,inp=inp)

        self._tTot += m.toc(t,"Couple ",self._is_time)

    def coupleProxi(self,mol,sio,orient,rate,inp="percent"):
        """Public method of :func:`_coupleProxi` for a adding specified rates of
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
        t = m.tic()
        if mol is not None: self._coupleProxi(mol,sio,orient,rate,inp=inp)
        self._tTot += m.toc(t,"Couple ",self._is_time)

    def finish(self,slx=0,isRand=True,isMol=False,isProps=True):
        """Finalize pore by adding siloxan bridges, filling empty bond with
        silanol groups, connecting geminal molecules, removing marked atoms,
        converting the grid to individual molecules and finally moving the pore
        into position. Furthermore the allocation is printed if needed.

        Parameters
        ----------
        slx : float
            Rate of siloxcan bridges
        isRand : bool
            True to randomly distribute rounding error charge on block oxygen atoms
        isMol : bool
            True to calculate allocation in micro molar
        isProps : bool
            True to calculate pore properties
        """
        # Start timer
        t = m.tic()

        # Analysis
        alloc = self._allocation(isMol=isMol) # Calculate allocation
        rough = self._rough()                 # Calculate pore roughness

        # Finalize
        self._siloxan(slx)      # Create silixan bridges
        self._silanolParallel() # Fill empty sites with silanol
        self._connect()         # Connect geminal molecules into one
        self._bonding.delete()  # Delete marked atoms
        self._objectify()       # Move all silica and oxygenes into unique mols
        self._sort()            # Sort molecules
        self._excess(isRand)    # Distribute excess charge
        self._position()        # Position the pore considering pbc
        self._overlap()         # Check for silica or oxygen atoms overlapping

        self._tTot += m.toc(t,"Finish ",self._is_time)

        if self._is_time:
            print("----------------------------")
            print("Total   - runtime = "+"%6.3f"%self._tTot+" s")
            print()

        if isProps:
            # Initialize
            diam  = self._diam
            size  = self._size
            site  = self._site
            drill = self._drill

            # Number of groups
            numPSL  = sum([1 for x in site if x[4]==0 and x[5] is None])
            numOSL  = sum([1 for x in site if x[4]==1 and x[5] is None])
            numPSLG = sum([1 for x in site if x[4]==0 and x[5] is not None])
            numOSLG = sum([1 for x in site if x[4]==1 and x[5] is not None])

            # Volume and surface
            vol   = math.pi*diam**2/4*size[2]
            surfP = math.pi*diam*size[2]
            surfO = size[0]*size[1]-math.pi*diam

            # Hydroxylation
            hydroP = (numPSL+numPSLG*2)/surfP
            hydroO = (numOSL+numOSLG*2)/surfO/2

            # Dimensions
            if   drill=="x": sizeStr = "%4.1f"%size[2]+" "+"%4.1f"%size[1]+" "+"%4.1f"%size[0]
            elif drill=="y": sizeStr = "%4.1f"%size[0]+" "+"%4.1f"%size[2]+" "+"%4.1f"%size[1]
            elif drill=="z": sizeStr = "%4.1f"%size[0]+" "+"%4.1f"%size[1]+" "+"%4.1f"%size[2]

            # Create dictionary
            props = {
                "D":   ["["+sizeStr+"]","nm"],
                "d":   ["%4.2f"%diam,   "nm"],
                "Rq":  ["%5.3f"%rough,  "nm"],
                "S":   ["%6.2f"%surfP,  "nm^2"],
                "V":   ["%6.2f"%vol,    "nm^3"],
                "hi":  ["%4.2f"%hydroP, "OH/nm^2"],
                "ho":  ["%4.2f"%hydroO, "OH/nm^2"],
                "hio": ["%5.3f"%(hydroP/hydroO),""],
                "Ni":  ["%4i"  %numPSL, "#"],
                "Nig": ["%4i"  %numPSLG,"#"],
                "Nii": ["%5.2f"%(numPSL/numPSLG),""],
                "No":  ["%4i"  %numOSL, "#"],
                "Nog": ["%4i"  %numOSLG,"#"],
                "Noo": ["%5.2f"%(numOSL/numOSLG),""],
                "ai":  ["%5.2f"%(alloc["pore"] ["rel"]*100),"%"],
                "ao":  ["%5.2f"%(alloc["drill"]["rel"]*100),"%"],
                "aio": ["%5.3f"%(alloc["pore"]["rel"]/alloc["drill"]["rel"]),""] if alloc["drill"]["rel"]>0 else ["%5.3f"%0,""],
                "C":   ["%7.5e"%self.getCharge(),"C"],
                "t":   ["%5.2f"%self._tTot,      "s"]
            }
            self._props = props

            # Print allocation
            unit = "[mumol/m^2]" if isMol else "[#/nm^2]"
            print("Surface allocation in "+unit)
            print("Unmodified Pore  - "+"%6.3f"% alloc["pore"]["tot"]     +",  Out - "+"%6.3f"% alloc["drill"]["tot"])
            print("Modified   Pore  - "+"%6.3f"% alloc["pore"]["mod"]     +",  Out - "+"%6.3f"% alloc["drill"]["mod"])
            print("Silanol    Pore  - "+"%6.3f"% alloc["pore"]["oh"]      +",  Out - "+"%6.3f"% alloc["drill"]["oh"])
            print("Relative   Pore  - "+"%6.3f"%(alloc["pore"]["rel"]*100)+"%, Out - "+"%6.3f"%(alloc["drill"]["rel"]*100)+"%")
            print()

            print("Properties")
            print("Dimensions            - "+props["D"] [0]+" "+props["D"]  [1])
            print("Diameter              - "+props["d"] [0]+" "+props["d"]  [1])
            print("Roughness             - "+props["Rq"][0]+" "+props["Rq"] [1])
            print("Surface               - "+props["S"] [0]+" "+props["S"]  [1])
            print("Volume                - "+props["V"] [0]+" "+props["V"]  [1])
            print("Hydroxylation in/out  - "+props["hi"][0]+"/"+props["ho"] [0]+" "+props["hi"][1]+" = "+props["hio"][0])
            print("Number of pore SL/SLG - "+props["Ni"][0]+"/"+props["Nig"][0]+" "+props["Ni"][1]+" = "+props["Nii"][0])
            print("Number of out  SL/SLG - "+props["No"][0]+"/"+props["Nog"][0]+" "+props["No"][1]+" = "+props["Noo"][0])
            print("Relative allocation   - "+props["ai"][0]+"/"+props["ao"] [0]+" "+props["ai"][1]+" = "+props["aio"][0])
            print("Excess charge         - "+props["C"] [0]+" "+props["C"]  [1])
            print("Total time            - "+props["t"] [0]+" "+props["t"]  [1])
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
            if not i==2: res *= size[i]

        pore = math.pi*(self._diam/2)**2

        return 2*res+pore


    ##################
    # Setter Methods #
    ##################
    def setGrid(self,grid):
        """Set the grid molecule list.

        Parameters
        ----------
        grid : list
            List of grid molecule identifier
        """
        self._grid = grid if isinstance(grid,list) else [grid]


    ##################
    # Getter Methods #
    ##################
    def getTime(self):
        """Return the total creation time.

        Returns
        -------
        time : float
            Used pore creation time
        """
        return self._tTot

    def getSize(self):
        """Return the pore size.

        Returns
        -------
        size : list
            Pore dimension
        """
        return self._size

    def getDiam(self):
        """Return the pore diameter.

        Returns
        -------
        diam : float
            Pore diameter
        """
        return self._diam

    def getReservoir(self):
        """Return the reservoir length.

        Returns
        -------
        res : float
            Reservoir length
        """
        return self._res

    def getGrid(self):
        """Return the grid molecules.

        Returns
        -------
        grid : list
            List of grid molecule identifier
        """
        return self._grid

    def getQ(self):
        """Return the adjusted grid charges.

        Returns
        -------
        q : dict
            Dictionary of grid molecules and their charges
        """
        return self._q

    def getDrill(self):
        """Return the drilling axis.

        Returns
        -------
        drill : str
            Drilling axis
        """
        return self._drill

    def getGap(self):
        """Return the gap vector.

        Returns
        -------
        gap : list
            Gap vector
        """
        return self._gap

    def getCenter(self):
        """Return the pore focal point.

        Returns
        -------
        center : list
            Pore focal point
        """
        return self._center

    def getProps(self):
        """Return the final pore properties.

        Returns
        -------
        props : dict
            Dictionary containing properties of the generated pore
        """
        return self._props
