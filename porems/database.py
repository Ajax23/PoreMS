################################################################################
# Database Class                                                               #
#                                                                              #
"""Very simple molecule database helper methods."""
################################################################################

# Create masses dictionary
masses = {"H":    1.0079,  # Hydrogen
          "He":   4.0026,  # Helium
          "Li":   6.9410,  # Lithium
          "Be":   9.0122,  # Beryllium
          "B":   10.8110,  # Boron
          "C":   12.0107,  # Carbon
          "N":   14.0067,  # Nitrogen
          "O":   15.9994,  # Oxygen
          "F":   18.9984,  # Fluorine
          "Ne":  20.1797,  # Neon
          "Na":  22.9897,  # Sodium
          "Mg":  24.3050,  # Magnesium
          "Al":  26.9815,  # Aluminum
          "Si":  28.0855,  # Silicon
          "P":   30.9738,  # Phosphorus
          "S":   32.0650,  # Sulfur
          "Cl":  35.4530,  # Chlorine
          "K":   39.0983,  # Potassium
          "Ar":  39.9480,  # Argon
          "Ca":  40.0780,  # Calcium
          "Sc":  44.9559,  # Scandium
          "Ti":  47.8670,  # Titanium
          "V":   50.9415,  # Vanadium
          "Cr":  51.9961,  # Chromium
          "Mn":  54.9380,  # Manganese
          "Fe":  55.8450,  # Iron
          "Ni":  58.6934,  # Nickel
          "Co":  58.9332,  # Cobalt
          "Cu":  63.5460,  # Copper
          "Zn":  65.3900,  # Zinc"
          "Ga":  69.7230,  # Gallium
          "Ge":  72.6400,  # Germanium
          "As":  74.9216,  # Arsenic
          "Se":  78.9600,  # Selenium
          "Br":  79.9040,  # Bromine
          "Kr":  83.8000,  # Krypton
          "Rb":  85.4678,  # Rubidium
          "Sr":  87.6200,  # Strontium
          "Y":   88.9059,  # Yttrium
          "Zr":  91.2240,  # Zirconium
          "Nb":  92.9064,  # Niobium
          "Mo":  95.9400,  # Molybdenum
          "Tc":  98.0000,  # Technetium
          "Ru": 101.0700,  # Ruthenium
          "Rh": 102.9055,  # Rhodium
          "Pd": 106.4200,  # Palladium
          "Ag": 107.8682,  # Silver
          "Cd": 112.4110,  # Cadmium
          "In": 114.8180,  # Indium
          "Sn": 118.7100,  # Tin
          "Sb": 121.7600,  # Antimony
          "I":  126.9045,  # Iodine
          "Te": 127.6000,  # Tellurium
          "Xe": 131.2930,  # Xenon
          "Cs": 132.9055,  # Cesium
          "Ba": 137.3270,  # Barium
          "La": 138.9055,  # Lanthanum
          "Ce": 140.1160,  # Cerium
          "Pr": 140.9077,  # Praseodymium
          "Nd": 144.2400,  # Neodymium
          "Pm": 145.0000,  # Promethium
          "Sm": 150.3600,  # Samarium
          "Eu": 151.9640,  # Europium
          "Gd": 157.2500,  # Gadolinium
          "Tb": 158.9253,  # Terbium
          "Dy": 162.5000,  # Dysprosium
          "Ho": 164.9303,  # Holmium
          "Er": 167.2590,  # Erbium
          "Tm": 168.9342,  # Thulium
          "Yb": 173.0400,  # Ytterbium
          "Lu": 174.9670,  # Lutetium
          "Hf": 178.4900,  # Hafnium
          "Ta": 180.9479,  # Tantalum
          "W":  183.8400,  # Tungsten
          "Re": 186.2070,  # Rhenium
          "Os": 190.2300,  # Osmium
          "Ir": 192.2170,  # Iridium
          "Pt": 195.0780,  # Platinum
          "Au": 196.9665,  # Gold
          "Hg": 200.5900,  # Mercury
          "Tl": 204.3833,  # Thallium
          "Pb": 207.2000,  # Lead
          "Bi": 208.9804,  # Bismuth
          "Po": 209.0000,  # Polonium
          "At": 210.0000,  # Astatine
          "Rn": 222.0000,  # Radon
          "Fr": 223.0000,  # Francium
          "Ra": 226.0000,  # Radium
          "Ac": 227.0000,  # Actinium
          "Pa": 231.0359,  # Protactinium
          "Th": 232.0381,  # Thorium
          "Np": 237.0000,  # Neptunium
          "U":  238.0289,  # Uranium
          "Am": 243.0000,  # Americium
          "Pu": 244.0000,  # Plutonium
          "Cm": 247.0000,  # Curium
          "Bk": 247.0000,  # Berkelium
          "Cf": 251.0000,  # Californium
          "Es": 252.0000,  # Einsteinium
          "Fm": 257.0000,  # Fermium
          "Md": 258.0000,  # Mendelevium
          "No": 259.0000,  # Nobelium
          "Lr": 262.0000,  # Lawrencium
          "Rf": 261.0000,  # Rutherfordium
          "Db": 262.0000,  # Dubnium
          "Sg": 266.0000,  # Seaborgium
          "Bh": 264.0000,  # Bohrium
          "Hs": 277.0000}  # Hassium


##################
# Getter Methods #
##################
def get_mass(symbol):
    """Return the mass of the specified atom.

        Parameters
        ----------
        symbol : str
            Atom symbol

        Returns
        -------
        mass : float
            Atom mass in :math:`\\frac g{mol}`
    """
    if symbol in masses:
        return masses[symbol]
    else:
        print("DB: Atom name not found.")
        return
