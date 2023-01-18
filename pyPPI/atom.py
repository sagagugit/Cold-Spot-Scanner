"""
Class representing atom
"""
ResiduesCodes = {
    'ALA': 'A',
    'ARG': 'R',
    'ASN': 'N',
    'ASP': 'D',
    'CYS': 'C',
    'GLU': 'E',
    'GLN': 'Q',
    'GLY': 'G',
    'HIS': 'H',
    'ILE': 'I',
    'LEU': 'L',
    'LYS': 'K',
    'MET': 'M',
    'PHE': 'F',
    'PRO': 'P',
    'SER': 'S',
    'THR': 'T',
    'TRP': 'W',
    'TYR': 'Y',
    'VAL': 'V',
    'GLX': 'Q',  # FAKE
    'ASX': 'N',  # FAKE
    'UNK': 'T'  # SUPER-FAKE
}


class atom(object):
    def __init__(self, atom_or_hetem, atom_num, atom_symbol, residue, chain, resId, iCode, atom_type, coord, occupancy, atom_index=0, beta_factor=0):
        self.neighbors = []
        self.atom_or_hetem = atom_or_hetem
        self.symbol = atom_symbol
        self.residue = residue
        self.chain = chain
        self.resId = int(resId)
        self.iCode = iCode
        self.x, self.y, self.z = coord
        self.coord = coord
        self.atomType = atom_type
        self.atomNum = atom_num
        self.occupancy = occupancy
        # used by pdb:
        self.atomIndex = atom_index
        self.tempFactor = beta_factor
        self.pseudoChain = None  # instead of chain

    def distance(self, other_atom):
        """Computes the distance to other atom
        :param other_atom: atom to compute distance to
        :return: distance
        """
        return (self.x - other_atom.x) ** 2 + (self.y - other_atom.y) ** 2 + (self.z - other_atom.z) ** 2

    def distanceFromXYZ(self, XYZ):
        """Computes distance from a coordinate XYZ
        :param XYZ:
        :return:
        """
        return (self.x - XYZ[0]) ** 2 + (self.y - XYZ[1]) ** 2 + (self.z - XYZ[2]) ** 2

    def resCode(self):
        global ResiduesCodes
        return ResiduesCodes[self.residue]

    def __str__(self):
        return self.chain + str(self.resId) + "|" + self.atomNum


class water(atom):
    """
    Water molecule
    """

    def __init__(self, atom_or_hetem, atom_num, atom_symbol, residue, chain, resId, iCode, atom_type, coord, occupancy):
        super(water, self).__init__(atom_or_hetem, atom_num, atom_symbol, residue, chain, resId, iCode, atom_type, coord, occupancy)
