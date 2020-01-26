################################################################################
# Essential Pack                                                               #
#                                                                              #
"""This file contains essential molecule classes and serves as an exemplary
molecule generation process."""
################################################################################


from porems.molecule import Molecule


class Alkane(Molecule):
    """Using this class linear alkane molecules can be easily constructed with
    the only input of the carbon atom number. Hydrogens are then automatically
    added.

    This class extends the molecule class :class:`porems.molecule.Molecule`.

    Parameters
    ----------
    length : integer
        Number of carbon atoms
    name : string, None, optional
        Molecule name
    short : string, None, optional
        Molecule short name
    is_h : bool, optional
        True if hydrogens are needed
    """
    def __init__(self, length, name=None, short=None, is_h=True):
        # Call super class
        super(Alkane, self).__init__()

        # Define molecule names
        if name is not None:
            self.set_name(name)
        if short is not None:
            self.set_short(short)

        # Define bond lengths and angles
        b = {"cc": 0.153, "ch": 0.109}
        a = {"ccc": 30.00, "cch": 109.47}

        # Add carbons
        self.add("C", [0, 0, 0])

        angle = a["ccc"]
        for i in range(length-1):
            angle *= -1
            self.add("C", self.get_num()-1, r=b["cc"], theta=angle)

        # Add hydrogens
        if is_h:
            if length > 1:
                # Boundaries
                end = self.get_num()-1

                for i in range(3):
                    self.add("H", 0, bond=[0, 1], r=b["ch"], theta=a["cch"], phi=60+120*i)
                    self.add("H", end, bond=[end, end-1], r=b["ch"], theta=a["cch"], phi=60+120*i)

                # Interior
                angle = -90
                for i in range(1, length-1):
                    angle *= -1
                    self.add("H", i, r=b["ch"], theta=angle, phi=a["cch"])
                    self.add("H", i, r=b["ch"], theta=angle, phi=-a["cch"])

            # Methane
            else:
                self.add("H", 0, r=b["ch"])
                for i in range(3):
                    self.add("H", 0, r=b["ch"], theta=a["cch"], phi=i*120)

        # Move to zero
        self.zero()


class Alcohol(Molecule):
    """Using this class linear alcohol molecules can be easily constructed with
    the only input of the carbon atom number. Hydrogens and the hydroxy group
    are then automatically added.

    This class extends the molecule class :class:`porems.molecule.Molecule`.

    Parameters
    ----------
    length : integer
        Number of carbon atoms
    name : string, None, optional
        Molecule name
    short : string, None, optional
        Molecule short name
    is_h : bool, optional
        True if hydrogens are needed
    """
    def __init__(self, length, name=None, short=None, is_h=True):
        # Call super class
        super(Alcohol, self).__init__()

        # Define molecule names
        if name is not None:
            self.set_name(name)
        if short is not None:
            self.set_short(short)

        # Define bond lengths and angles
        b = {"cc": 0.153, "ch": 0.109, "co": 0.143, "oh": 0.098}
        a = {"ccc": 30.00, "cch": 109.47, "occ": 30.00, "coh": 109.47}

        # Add carbons
        self.add("C", [0, 0, 0])

        angle = a["ccc"]
        for i in range(length-1):
            angle *= -1
            self.add("C", self.get_num()-1, r=b["cc"], theta=angle)

        # Add hydroxy
        self.add("O", self.get_num()-1, r=b["co"], theta=-angle)
        self.add("H", self.get_num()-1, r=b["oh"], theta=angle)

        # Add hydrogens
        if is_h:
            if length > 1:
                # Boundary
                for i in range(3):
                    self.add("H", 0, bond=[0, 1], r=b["ch"], theta=a["cch"], phi=120*i)

                # Interior
                angle = -90
                for i in range(1, length):
                    angle *= -1
                    self.add("H", i, r=b["ch"], theta=angle, phi=a["cch"])
                    self.add("H", i, r=b["ch"], theta=angle, phi=-a["cch"])

            # Methanol
            else:
                for i in range(3):
                    self.add("H", 0, r=b["ch"], theta=a["cch"], phi=i*120)

        # Move to zero
        self.zero()


class Ketone(Molecule):
    """Using this class linear ketone molecules can be easily constructed with
    the only input of the carbon atom number and oxygen position.
    Hydrogens are then automatically added.

    This class extends the molecule class :class:`porems.molecule.Molecule`.

    Parameters
    ----------
    length : integer
        Number of carbon atoms
    pos : integer
        Position of the oxygen atom
    name : string, None, optional
        Molecule name
    short : string, None, optional
        Molecule short name
    is_h : bool
        True if hydrogens are needed
    """
    def __init__(self, length, pos, name=None, short=None, is_h=True):
        # Call super class
        super(Ketone, self).__init__()

        # Check input
        if length < 3:
            print("Specified length is too small for ketones ...")
            return

        # Define molecule names
        if name is not None:
            self.set_name(name)
        if short is not None:
            self.set_short(short)

        # Define bond lengths and angles
        b = {"cc": 0.153, "ch": 0.109, "co": 0.123}
        a = {"ccc": 30.00, "cch": 109.47}

        # Add carbons
        self.add("C", [0, 0, 0])

        angle = a["ccc"]
        for i in range(length-1):
            angle *= -1
            self.add("C", self.get_num()-1, r=b["cc"], theta=angle)

        # Add oxygen
        angle = -90 if pos % 2 == 0 else 90
        self.add("O", pos-1, r=b["co"], theta=angle)

        # Add hydrogens
        if is_h:
            # Boundaries
            end = self.get_num()-1

            for i in range(3):
                self.add("H", 0, bond=[0,  1], r=b["ch"], theta=a["cch"], phi=60+120*i)
                self.add("H", end-1, bond=[end-1, end-2], r=b["ch"], theta=a["cch"], phi=60+120*i)

            # Interior
            angle = -90
            for i in range(1, length-1):
                angle *= -1
                if not i == pos-1:
                    self.add("H", i, r=b["ch"], theta=angle, phi=a["cch"])
                    self.add("H", i, r=b["ch"], theta=angle, phi=-a["cch"])

        # Move to zero
        self.zero()


class TMS(Molecule):
    """This class defines trimethylsilyl (TMS) bound to a silicon-grid.

    This class extends the molecule class :class:`porems.molecule.Molecule`.

    Parameters
    ----------
    name : string, optional
        Molecule name
    short : string, optional
        Molecule short name
    charge : float, optional
        Excess charge of the entire grid molecule
    compress : float, optional
        Compress molecule sidechain angles
    is_si : bool, optional
        True if the terminus should be a lone silicon, False for a CH3 group
    """
    def __init__(self, name="tms", short="TMS", charge=1.070428, compress=0, is_si=True):
        # Call super class
        super(TMS, self).__init__()

        # Define molecule names
        self.set_name(name)
        self.set_short(short)

        # Check silicon
        si = "Si" if is_si else "Ci"

        # Define bond lengths and angles
        b = {"sio": 0.155, "sic": 0.186, "ch": 0.109}
        a = {"ccc": 30.00, "cch": 109.47, "siosi": 126.12}

        # Build silyl chain
        self.add(si, [0, 0, 0])
        self.add("O", 0, r=b["sio"])
        self.add(si, 1, bond=[1, 0], r=b["sio"], theta=a["siosi"])

        # Add methyl
        for i in range(3):
            self.add("C", 2, bond=[2, 1], r=b["sic"], theta=a["cch"]+compress, phi=60+120*i)

        # Add hydrogens
        for i in range(3, 5+1):
            for j in range(3):
                self.add("H", i, bond=[i, 2], r=b["ch"], theta=a["cch"]+compress, phi=60+120*j)

        # If not silicon ending
        for i in range(3):
            if not is_si:
                self.add("H", 0, bond=[0, 1], r=b["ch"], theta=a["cch"], phi=60+120*i)

        # Set charge
        self.set_charge(charge)
