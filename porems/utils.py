################################################################################
# Utils                                                                        #
#                                                                              #
"""Here popular basic methods are noted."""
################################################################################


import os
import time
import pickle
import fileinput

from shutil import copyfile


def mkdirp(directory):
    """Create directory if it does not exist.

    Parameters
    ----------
    directory : string
        Directory name
    """
    if not os.path.exists(directory):
        os.makedirs(directory)


def copy(source, target):
    """Copy a specified file to a specified location.

    Parameters
    ----------
    source : string
        Link to requested file
    target : string
        Link to requested new file
    """
    copyfile(source, target)


def column(data):
    """Convert given row list matrix into column list matrix

    Parameters
    ----------
    data : list
        Row data matrix

    Returns
    -------
    data_col : list
        column data matrix
    """
    num_row = len(data)
    num_col = len(data[0])

    data_col = [[] for i in range(num_col)]

    for i in range(num_row):
        for j in range(num_col):
            data_col[j].append(data[i][j])

    return data_col


def tic():
    """MATLAB tic replica - return current time.

    Returns
    -------
    time : float
        Current time in seconds
    """
    return time.time()


def toc(t, message="", is_print=True):
    """MATLAB toc replica - return time difference to tic and alternatively
    print a message.

    Parameters
    ----------
    t : float
        Starting time - given by tic function
    message : string, optional
        Custom output message
    is_print : bool, optional
        True for printing an output message

    Returns
    -------
    time : float
        Time difference
    """
    if message:
        message += " - runtime = "

    t_diff = time.time()-t

    if is_print:
        print(message+"%6.3f" % t_diff+" s")

    return t_diff


def replace(file_link, old, new):
    """Replace a given string in a file.

    Parameters
    ----------
    file_link : string
        Link to requested file
    old : string
        String to be replaced in file
    new : string
        New string to be written
    """
    for line in fileinput.input(file_link, inplace=True):
        print(line.rstrip().replace(old, new))


def save(obj, link):
    """Save an object using pickle in the specified link.

    Parameters
    ----------
    obj : Object
        Object to be saved
    link : string
        Specific link to save object
    """
    with open(link, "wb") as f:
        pickle.dump(obj, f)


def load(link):
    """Load pickled object from the specified folder.

    Parameters
    ----------
    link : string
        Specific link to load object

    Returns
    -------
    obj : Object
        Loaded object
    """
    with open(link, 'rb') as f:
        return pickle.load(f)


def mumol_m2_to_mols(c, A):
    """Convert the concentration in :math:`\\frac{\\mu\\text{mol}}{\\text{m}^2}`
    to number of molecules.

    The concentration is given as

    .. math::

        c=\\frac{N}{N_A\\cdot A}

    with surface :math:`A` and avogadro constant
    :math:`N_A=6.022\\cdot10^{23}\\,\\text{mol}^{-1}`. Assuming that the unit of
    the concentration is :math:`\\mu\\text{mol}\\cdot\\text{m}^{-2}` and of the
    surface is :math:`\\text{nm}^2`, the number of molecules is given by

    .. math::

        N=c\\cdot N_A\\cdot A
        =[c]\\cdot10^{-24}\\,\\frac{\\text{mol}}{\\text{nm}^2}\\cdot6.022\\cdot10^{23}\\,\\frac{1}{\\text{mol}}\\cdot[A]\\,\\text{nm}^2

    and thus

    .. math::

        \\boxed{N=0.6022\\cdot[c]\\cdot[A]}\\ .

    Parameters
    ----------
    c : float
        Concentration in :math:`\\frac{\\mu\\text{mol}}{\\text{m}^2}`
    A : float
        Surface in :math:`\\text{nm}^2`

    Returns
    -------
    N : float
        Number of molecules
    """
    return 0.6022*c*A


def mols_to_mumol_m2(N, A):
    """Convert the number of molecules to concentration in
    :math:`\\frac{\\mu\\text{mol}}{\\text{m}^2}`.

    The concentration is given as

    .. math::

        c=\\frac{N}{N_A\\cdot A}

    with surface :math:`A` and avogadro constant
    :math:`N_A=6.022\\cdot10^{23}\\,\\text{mol}^{-1}`. Assuming that the unit of
    the concentration is :math:`\\mu\\text{mol}\\cdot\\text{m}^{-2}` and of the
    surface is :math:`\\text{nm}^2`, the number of molecules is given by

    .. math::

        c=\\frac{[N]}{6.022\\cdot10^{23}\\,\\text{mol}^{-1}\\cdot[A]\\,\\text{nm}^2}

    and thus

    .. math::

        \\boxed{c=\\frac{N}{0.6022\\cdot[A]}\\cdot\\frac{\\mu\\text{mol}}{\\text{m}^2}}\\ .

    Parameters
    ----------
    N : int
        Number of molecules
    A : float
        Surface in :math:`\\text{nm}^2`

    Returns
    -------
    c : float
        Concentration in :math:`\\frac{\\mu\\text{mol}}{\\text{m}^2}`
    """
    return N/0.6022/A


def mmol_g_to_mumol_m2(c, SBET):
    """Convert the concentration :math:`\\frac{\\text{mmol}}{\\text{g}}`
    to :math:`\\frac{\\mu\\text{mol}}{\\text{m}^2}`.

    This is done by dividing the concentration per gram by the material surface
    per gram property :math:`S_\\text{BET}`

    .. math::

        \\boxed{c_A=\\frac{c_g}{S_\\text{BET}}\\cdot10^3}\\ .

    Parameters
    ----------
    c : float
        Concentration in :math:`\\frac{\\text{mmol}}{\\text{g}}`
    SBET : float
        Material surface in :math:`\\frac{\\text{m}^2}{\\text{g}}`

    Returns
    -------
    c : float
        Concentration in :math:`\\frac{\\mu\\text{mol}}{\\text{m}^2}`
    """
    return c/SBET*1e3


def mmol_l_to_mols(c, V):
    """Convert the concentration in :math:`\\frac{\\text{mmol}}{\\text{l}}`
    to number of molecules.

    The concentration in regard to volume is calculated by

    .. math::

        c=\\frac{N}{N_A\\cdot V}

    with volume :math:`V`. Assuming that the unit of the concentration is
    :math:`\\text{mmol}\\cdot\\text{l}` and of the volume is
    :math:`\\text{nm}^3`, the number of molecules is given by

    .. math::

        N=c\\cdot N_A\\cdot V
        =[c]\\cdot10^{-27}\\,\\frac{\\text{mol}}{\\text{nm}^3}\\cdot6.022\\cdot10^{23}\\,\\frac{1}{\\text{mol}}\\cdot[V]\\,\\text{nm}^3

    and thus

    .. math::

        \\boxed{N=6.022\\cdot10^{-4}\\cdot[c]\\cdot[V]}\\ .

    Parameters
    ----------
    c : float
        Concentration in :math:`\\frac{\\text{mmol}}{\\text{l}}`
    V : float
        Surface in :math:`\\text{nm}^3`

    Returns
    -------
    N : float
        Number of molecules
    """
    return 6.022e-4*c*V


def mols_to_mmol_l(N, V):
    """Convert the number of molecules to concentration in
    :math:`\\frac{\\text{mmol}}{\\text{l}}`.

    The concentration in regard to volume is calculated by

    .. math::

        c=\\frac{N}{N_A\\cdot V}

    with volume :math:`V`. Assuming that the unit of the concentration is
    :math:`\\text{mmol}\\cdot\\text{l}` and of the volume is
    :math:`\\text{nm}^3`, the number of molecules is given by

    .. math::

        c=\\frac{N}{6.022\\cdot10^{23}\\cdot[V]\\,\\text{nm}^3}

    and thus

    .. math::

        \\boxed{c=\\frac{N}{6.022\\times10^{-4}\\cdot[V]}\\cdot\\frac{\\text{mmol}}{\\text{l}}}\\ .

    Parameters
    ----------
    N : int
        Number of molecules
    V : float
        Surface in :math:`\\text{nm}^3`

    Returns
    -------
    c : float
        Concentration in :math:`\\frac{\\text{mmol}}{\\text{l}}`
    """
    return N/6.022e-4/V
