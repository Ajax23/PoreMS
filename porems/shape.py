################################################################################
# Shape Pack                                                                   #
#                                                                              #
"""This file contains shape definitions to be cut out from the crystal block.
"""
################################################################################


import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

import porems.geometry as geometry


class Shape():
    """This class is a container for individual shape classes.

    Parameters
    ----------
    inp : dictionary
        Dictionary of necessary inputs
    """
    def __init__(self, inp):
        self._inp = inp

        # Calculate angle and normal vector for rotation
        self._angle = geometry.angle(inp["central"], geometry.main_axis("z"))
        self._normal = geometry.cross_product(inp["central"], geometry.main_axis("z"))

        # Calculate distance towards central axis start
        self._dist_start = geometry.vector(self._centroid, inp["centroid"])
        self._dist_zero = geometry.vector(geometry.rotate(inp["centroid"], self._normal, -self._angle, True), self._centroid)


    ##################
    # Helper Methods #
    ##################
    def convert(self, data, to_zero=True):
        """This helper method converts the given data to match the global
        cenrtal axis and its start.

        Parameters
        ----------
        data : list
            Data input
        to_zero : bool
            True to convert data towards zero axis, False to convert data from
            zero axis to central axis.

        Returns
        -------
        data : list
            Converted input
        """
        # Rotate towards main axis ot to zero axis
        data = geometry.rotate(data, self._normal, -self._angle if to_zero else self._angle, True)

        # Translate to zero or to start
        dist = self._dist_zero if to_zero else self._dist_start
        data = [data[i]+dist[i] for i in range(3)]

        return data

    def plot(self, z=0, num=100):
        """Plot surface and rim.

        Parameters
        ----------
        z : float, optional
            Position on the axis
        num : integer, optional
            Number of points
        """
        fig = plt.figure()

        # Surface
        ax = fig.add_subplot(121, projection="3d")
        ax.plot_surface(*self.surf(num=100))

        # Rim
        ax = fig.add_subplot(122, projection="3d")
        ax.plot3D(*[x[0] for x in self.rim(z, num)])

        # Normal
        # import porems.utils as utils
        # vec = [3.17290646, 4.50630614, 0.22183271]
        # line = [self.convert([0, 0, 0], False), vec,
        #         self.convert(self.normal(vec), False),
        #         self.convert([x*6 for x in geometry.unit(self.normal(vec))], False)]
        # ax.plot3D(*utils.column(line))

        plt.show()


class Cylinder(Shape):
    """This class defines a cylindric shape. Needed inputs are

    * **central** - Central axis
    * **centroid** - Centroid of block
    * **len_block** - Length of block
    * **len_cyl** - Cylinder length
    * **diameter** - Cylinder diameter

    Parameters
    ----------
    inp : dictionary
        Dictionary of necessary inputs
    """
    def __init__(self, inp):
        # Set centroid
        self._centroid = [0, 0, inp["len_cyl"]/2]

        # Call super class
        super(Cylinder, self).__init__(inp)


    ############
    # Function #
    ############
    def Phi(self, r, phi, z):
        """Surface function.

        Parameters
        ----------
        r : float
            Radius
        phi : float
            Polar angle
        z : float
            Distance ins z-axis

        Returns
        -------
        pos : list
            Cartesian coordinates for given polar angle and z-distance
        """
        # Calculate coordinates
        x = np.outer(r, np.cos(phi))
        y = np.outer(r, np.sin(phi))
        z = z if isinstance(z, list) else [z]
        z = np.outer(z, np.ones(len(z)))

        return self.convert([x, y, z], False)

    def d_Phi_phi(self, r, phi, z):
        """Derivative of the surface function considering the polar angle.

        Parameters
        ----------
        r : float
            Radius
        phi : float
            Polar angle
        z : float
            Distance ins z-axis

        Returns
        -------
        pos : list
            Cartesian coordinates for given polar angle and z-distance
        """
        x = -r*np.sin(phi)
        y = r*np.cos(phi)
        z = 0

        return [x, y, z]

    def d_Phi_z(self, r, phi, z):
        """Derivative of the surface function considering the z-axis.

        Parameters
        ----------
        r : float
            Radius
        phi : float
            Polar angle
        z : float
            Distance ins z-axis

        Returns
        -------
        pos : list
            Cartesian coordinates for given polar angle and z-distance
        """
        x = 0
        y = 0
        z = 1

        return [x, y, z]


    ############
    # Features #
    ############
    def normal(self, pos):
        """Calculate unit normal vector on surface for given position.

        Parameters
        ----------
        pos : list
            Position

        Returns
        -------
        normal : list
            Unit normal vector
        """
        # Initialize
        x, y, z = self.convert(pos)

        # Cartesian to polar
        r = math.sqrt(x**2+y**2)
        phi = geometry.angle_polar([x, y, z])
        z = z

        # Calculate derivatives
        d_Phi_phi = self.d_Phi_phi(r, phi, z)
        d_Phi_z = self.d_Phi_z(r, phi, z)

        # Calculate normal vector
        return geometry.cross_product(d_Phi_phi, d_Phi_z)

    def is_in(self, pos):
        """Check if given position is inside of shape. Needed inputs are

        Parameters
        ----------
        pos : list
            Position

        Returns
        -------
        is_in : bool
            True if position is inside of shape
        """
        # Check if within shape
        if geometry.length(self.normal(pos)) < self._inp["diameter"]/2:
            pos_zero = self.convert(pos)
            return pos_zero[2]>0 and pos_zero[2]<self._inp["len_cyl"]
        else:
            return False


    #########
    # Shape #
    #########
    def rim(self, z, num=100):
        """Return x and y values for given z-position.

        Parameters
        ----------
        z : float
            Position on the axis
        num : integer, optional
            Number of points

        Returns
        -------
        positions : list
            x and y arrays of the surface rim on the z-position
        """
        phi = np.linspace(0, 2*np.pi, num)
        r = self._inp["diameter"]/2

        return self.Phi(r, phi, z)

    def surf(self, num=100):
        """Return x, y and z values for the shape.

        Parameters
        ----------
        num : integer, optional
            Number of points

        Returns
        -------
        positions : list
            x, y and z arrays of the surface rim
        """
        phi = np.linspace(0, 2*np.pi, num)
        r = np.ones(num)*self._inp["diameter"]/2
        z = np.linspace(0, self._inp["len_cyl"], num)

        return self.Phi(r, phi, z)


    ##############
    # Properties #
    ##############
    def volume(self):
        """Calculate volume.

        Returns
        -------
        volume : float
            Volume
        """
        return math.pi*(self._inp["len_cyl"]/2)**2*self._inp["len_cyl"]

    def surface(self):
        """Calculate inner surface.

        Returns
        -------
        surface : float
            Inner surface
        """
        return 2*math.pi*self._inp["len_cyl"]/2*self._inp["len_cyl"]
