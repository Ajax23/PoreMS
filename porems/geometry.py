################################################################################
# Geometry                                                                     #
#                                                                              #
"""Here basic geometric functions are noted."""
################################################################################


import math


def dot_product(vec_a, vec_b):
    """Calculate the dot product of two vectors
    :math:`\\boldsymbol{a},\\boldsymbol{b}\\in\\mathbb{R}^n`

    .. math::

        \\text{dot}(\\boldsymbol{a},\\boldsymbol{b})=
        \\begin{pmatrix}a_1\\\\\\vdots\\\\a_n\\end{pmatrix}\\cdot
        \\begin{pmatrix}b_1\\\\\\vdots\\\\b_n\\end{pmatrix}=
        a_1\\cdot b_1+a_2\\cdot b_2+\\dots+a_n\\cdot b_n.

    Parameters
    ----------
    vec_a : list
        First vector :math:`\\boldsymbol{a}`
    vec_b : list
        Second vector :math:`\\boldsymbol{b}`

    Returns
    -------
    dot : float
        Dot product value
    """
    return sum((a*b) for a, b in zip(vec_a, vec_b))


def length(vec):
    """Calculate the length of a vector
    :math:`\\boldsymbol{a}\\in\\mathbb{R}^n`

    .. math::

        \\text{length}(\\boldsymbol{a})=|\\boldsymbol{a}|
        =\\sqrt{\\boldsymbol{a}\cdot\\boldsymbol{a}}

    Parameters
    ----------
    vec : list
        Vector a

    Returns
    -------
    length : float
        Vector length
    """
    return math.sqrt(dot_product(vec, vec))


def vector(pos_a, pos_b):
    """Calculate the vector between to two positions
    :math:`\\boldsymbol{a},\\boldsymbol{b}\\in\\mathbb{R}^n`

    .. math::

        \\text{vec}(\\boldsymbol{a},\\boldsymbol{b})
        =\\begin{pmatrix}b_1-a_1\\\\\\vdots\\\\b_n-a_n\\end{pmatrix}

    Parameters
    ----------
    pos_a : list
        First position :math:`\\boldsymbol{a}`
    pos_b : list
        Second position :math:`\\boldsymbol{b}`

    Returns
    -------
    vector : list
        Bond vector
    """
    # Check dimensions
    if not len(pos_a) == len(pos_b):
        print("Vector: Wrong dimensions...")
        return

    # Calculate vector
    return [pos_b[i]-pos_a[i] for i in range(len(pos_a))]


def unit(vec):
    """Transform a vector :math:`\\boldsymbol{a}\\in\\mathbb{R}^n` into a
    unit vector

    .. math::

        \\text{unit}(\\boldsymbol{a})
        =\\frac{\\boldsymbol{a}}{|\\boldsymbol{a}|}

    Parameters
    ----------
    vec : list
        Vector a

    Returns
    -------
    vec : list
        Vector
    """
    vec_length = length(vec)

    return [x/vec_length if not vec_length == 0 else x for x in vec]


def cross_product(vec_a, vec_b):
    """Calculate the cross product of two three-dimensional vectors
    :math:`\\boldsymbol{a},\\boldsymbol{b}\\in\\mathbb{R}^3`

    .. math::

        \\text{cross}(\\boldsymbol{a},\\boldsymbol{b})=\\begin{pmatrix}
        a_2\\cdot b_3-a_3\\cdot b_2\\\\
        a_3\\cdot b_1-a_1\\cdot b_4\\\\
        a_1\\cdot b_2-a_2\\cdot b_1
        \\end{pmatrix}

    Parameters
    ----------
    vec_a : list
        First vector :math:`\\boldsymbol{a}`
    vec_b : list
        Second vector :math:`\\boldsymbol{b}`

    Returns
    -------
    vec : list
        Cross product vector
    """
    vec = []
    vec.append(vec_a[1]*vec_b[2]-vec_a[2]*vec_b[1])
    vec.append(vec_a[2]*vec_b[0]-vec_a[0]*vec_b[2])
    vec.append(vec_a[0]*vec_b[1]-vec_a[1]*vec_b[0])

    return vec


def angle(vec_a, vec_b, is_deg=True):
    """Calculate the angle between two vectors
    :math:`\\boldsymbol{a},\\boldsymbol{b}\\in\\mathbb{R}^n`

    .. math::

        \\text{angle}=\\cos^{-1}\\frac{\\boldsymbol{a}\cdot\\boldsymbol{b}}
        {|\\boldsymbol{a}||\\boldsymbol{a}|}

    Parameters
    ----------
    vec_a : list
        First vector :math:`\\boldsymbol{a}`
    vec_b : list
        Second vector :math:`\\boldsymbol{b}`
    is_deg : bool, optional
        True if the output should be in degree

    Returns
    -------
    angle : float
        Angle
    """
    angle = math.acos(dot_product(vec_a, vec_b)/(length(vec_a)*length(vec_b)))

    return angle*180/math.pi if is_deg else angle


def angle_polar(pos, is_deg=False):
    """Calculate the polar angle of a position vector
    :math:`\\boldsymbol{a}\\in\\mathbb{R}^3`, which is the angle of the
    x-axis towards the reflected position vector on the x-y-plane

    .. math::

        \\text{polar}(\\boldsymbol{a})=\\arctan2(x,y)\\left\\{
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
        Position vector :math:`\\boldsymbol{a}`
    is_deg : bool, optional
        True if the output should be in degree

    Returns
    -------
    angle : float
        Polar angle
    """
    angle = math.atan2(pos[1], pos[0])

    return angle*180/math.pi if is_deg else angle


def angle_azi(pos, is_deg=False):
    """Calculate the azimuthal angle of a position vector
    :math:`\\boldsymbol{a}\\in\\mathbb{R}^3`, which is the angle of the
    position vector towards the x-y-plane

    .. math::

        \\text{azimut}(\\boldsymbol{a})
        =\\cos^{-1}\\frac{y}{|\\boldsymbol{a}|}

    with :math:`y` as the second vector entry.

    Parameters
    ----------
    pos : list
        Position vector :math:`\\boldsymbol{a}`
    is_deg : bool, optional
        True if the output should be in degree

    Returns
    -------
    angle : float
        Azimuthal angle
    """
    try:
        angle = math.acos(pos[2]/length(pos))
    except(ZeroDivisionError):
        angle = math.acos(0)

    return angle*180/math.pi if is_deg else angle


def main_axis(inp, dim=3):
    """Return the three-dimensional unit-vector of the main axes.
    Input is either integer or string

    * 1 or "x" for the x-axis
    * 2 or "y" for the y-axis
    * 3 or "z" for the z-axis

    Parameters
    ----------
    inp : integer, string
        Axis type input
    dim : integer, optional
        Number of dimensions

    Returns
    -------
    vec : list
        Unit vector
    """
    # Error message
    axis_error = "Wrong axis definition..."

    # Process input
    if isinstance(inp, str):
        if inp == "x":
            axis = 1
        elif inp == "y":
            axis = 2
        elif inp == "z":
            axis = 3
        else:
            return axis_error
    elif isinstance(inp, int):
        if inp == 1 or inp == 2 or inp == 3:
            axis = inp
        else:
            return axis_error
    else:
        return axis_error

    # Return vector
    return [1 if i == axis-1 else 0 for i in range(dim)]


def rotate(data, axis, angle, is_deg, dim=3):
    """Rotate a vector :math:`\\boldsymbol{a}\\in\\mathbb{R}^3`
    along an axis :math:`\\boldsymbol{b}\\in\\mathbb{R}^3` with angle
    :math:`\\alpha\\in\\mathbb{R}`. The rotation is performed using the
    rotation-matrix

    .. math::

        \\boldsymbol{R}_\\boldsymbol{n}(\\alpha)=\\begin{pmatrix}
        n_1^2 (1-\\cos\\alpha)+   \\cos\\alpha&n_1n_2(1-\\cos\\alpha)-n_3\\sin\\alpha&n_1n_3(1-\\cos\\alpha)+n_2\\sin\\alpha\\\\
        n_2n_1(1-\\cos\\alpha)+n_3\\sin\\alpha&n_2^2 (1-\\cos\\alpha)+   \\cos\\alpha&n_2n_3(1-\\cos\\alpha)-n_1\\sin\\alpha\\\\
        n_3n_1(1-\\cos\\alpha)-n_2\\sin\\alpha&n_3n_2(1-\\cos\\alpha)+n_1\\sin\\alpha&n_3^2 (1-\\cos\\alpha)+   \\cos\\alpha
        \\end{pmatrix}


    where :math:`n_i` are the entries for the unit vector
    :math:`\\boldsymbol{n}` of the axis. The new coordinates
    :math:`\\boldsymbol{c}` are then calculated using a matrix-vector
    multiplication

    .. math::

        \\boldsymbol{c}=\\boldsymbol{R}_\\boldsymbol{n}\\boldsymbol{a}.

    Parameters
    ----------
    data : list
        Vector :math:`\\boldsymbol{a}`
    axis : integer, string, list
        Axis :math:`\\boldsymbol{b}`
    angle : float
        Angle
    is_deg : bool
        True if the input is in degree
    dim : integer, optional
        Number of dimensions

    Returns
    -------
    coord : list
        Vector c as the result of the rotation
    """
    # Angle
    angle = angle*math.pi/180 if is_deg else angle

    # Set vector
    if isinstance(axis, list):
        if len(axis) == dim:
            n = axis
        elif len(axis) == 2:
            n = vector(axis[0], axis[1])
        else:
            print("Rotate: Wrong vector dimensions.")
            return
    else:
        n = main_axis(axis)
        if isinstance(n, str):
            print("Rotate: "+n)
            return

    # Unit vector
    n = unit(n)

    # Define rotation matrix
    n1 = n[0]
    n2 = n[1]
    n3 = n[2]

    c = math.cos(angle)
    s = math.sin(angle)

    r = [[n1*n1*(1.-c) + c, n1*n2*(1.-c) - n3*s,  n1*n3*(1.-c) + n2*s],
         [n2*n1*(1.-c) + n3*s, n2*n2*(1.-c) + c,  n2*n3*(1.-c) - n1*s],
         [n3*n1*(1.-c) - n2*s, n3*n2*(1.-c) + n1*s,  n3*n3*(1.-c) + c]]

    # Rotate
    return [data[0]*r[i][0]+data[1]*r[i][1]+data[2]*r[i][2] for i in range(dim)]
