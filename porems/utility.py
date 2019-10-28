################################################################################
# Utility                                                                      #
#                                                                              #
"""Here popular basic methods are noted. Furthermore terminal commands are
specified depending on the operating system."""
################################################################################


import os
import time
import pickle
import fileinput

from shutil import copyfile, copytree


def isWin():
    """Check operating system type.

    Returns
    -------
    system : bool
        True if operating system is windows
    """
    return os.name == "nt"

def mkdirp(directory):
    """Create directory if it does not exist.

    Parameters
    ----------
    directory : str
        Directory name
    """
    if not os.path.exists(directory):
        os.makedirs(directory)

def read(fileLink):
    """Read given file and return a data matrix.

    Parameters
    ----------
    fileLink : str
        Link to requested file

    Returns
    -------
    data : list
        data matrix extracted from file
    """
    # Initialize
    data  = []

    # Read input daat from file
    fileIn = open(fileLink,"r")
    for line in fileIn:
        data.append(line.split())
    fileIn.close()

    return data

def column(data):
    """Convert given row list matrix into column list matrix

    Parameters
    ----------
    data : list
        Row data matrix

    Returns
    -------
    dat : list
        column data matrix
    """
    numRow = len(data)
    numCol = len(data[0])

    dat = [[] for i in range(numCol)]

    for i in range(numRow):
        for j in range(numCol):
            dat[j].append(data[i][j])

    return dat

def getFolder(folderDir="."):
    """Get a list of folders in given directory.

    Parameters
    ----------
    folderDir : str
        Link to requested folder

    Returns
    -------
    folders : list
        List of folders
    """
    return [f for f in os.listdir(folderDir) if os.path.isdir(os.path.join(folderDir, f))]

def getFile(folderDir=".",dont=None):
    """Get a list of files in given directory.

    Parameters
    ----------
    folderDir : str
        Link to requested folder
    dont : list, str or None
        Filenames not to include in output

    Returns
    -------
    folders : list
        List of folders
    """
    if not isinstance(dont,list): dont = [dont]

    fileList = [f for f in os.listdir(folderDir) if f not in dont]

    return fileList

def copy(source,target):
    """Copy a specified file to a specified location.

    Parameters
    ----------
    source : str
        Link to requested file
    target : str
        Link to requested new file
    """
    copyfile(source,target)

def copyDir(source,target):
    """Copy a specified folder to a specified location if it does not exist.

    Parameters
    ----------
    source : str
        Link to requested folder
    target : str
        Link to requested new folder
    """
    if not os.path.isdir(target):
        copytree(source,target)

def hasNumber(inputString):
    """Check if string contains numbers.

    Parameters
    ----------
    inputString : str
        Given string

    Returns
    -------
    check : bool
        True if string contains number
    """
    return any(char.isdigit() for char in inputString)

def isInt(inputString):
    """Check if the input can be converted to an integer.

    Parameters
    ----------
    inputString : str
        Given string

    Returns
    -------
    check : bool
        True if string is an integer
    """
    try:
        int(inputString)
        return True
    except ValueError:
        return False

def tic():
    """Matlab tic reproduction - return current time.

    Returns
    -------
    time : float
        Current time in seconds
    """
    return time.time()

def toc(t,message=None,isPrint=True):
    """Matlab toc reproduction - return time difference to tic and alternatively
    print a message.

    Parameters
    ----------
    t : float
        Starting time - given by tic function
    message : str
        Custom output message
    isPrint : bool
        True for printing an output message

    Returns
    -------
    time : float
        Time difference
    """
    if message==None:
        message  = ""
    else:
        message += " - runtime = "

    tDiff = time.time()-t

    if isPrint: print(message+"%6.3f"%tDiff+" s")

    return tDiff

def replace(fileLink,old,new):
    """Replace a given string in a file.

    Parameters
    ----------
    fileLink : str
        Link to requested file
    old : str
        String to be replaced in file
    new : str
        New string to be written
    """
    for line in fileinput.input(fileLink,inplace=True):
          print(line.rstrip().replace(old,new))

def save(obj,link):
    """Save an object using pickle in the specified link.

    Parameters
    ----------
    obj : Object
        Object to be saved
    link : str
        Specific link to save object
    """
    with open(link, "wb") as f: pickle.dump(obj,f)

def load(link):
    """Load pickled object from the specified folder.

    Parameters
    ----------
    link : str
        Specific link to load object

    Returns
    -------
    obj : Object
        Loaded object
    """
    return pickle.load(open(link, "rb"))

def removeNum(inputString):
    """Remove numbers and then spaces from given string.

    Parameters
    ----------
    inputString : str
        Given string

    Returns
    -------
    output : str
        Output string without numbers
    """
    return "".join([i for i in inputString if not i.isdigit()]).strip()
