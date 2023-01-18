"""
Accessible surface area calculations
"""
#from __future__ import print_function

import logging
from math import cos, sin, sqrt, pi

from .pdbReader import PDBReader
from .resReduce import mapUnbound

R_WATER = 1.4
# signtificant distance
KNOWN_RADIUS = {
    'C': 1.700,
    'N': 1.550,
    'O': 1.520,
    'S': 1.800,
    'H': 1.200,
    'X': 1.700  # for undetermined atom assume Carbon
}
SPIRAL_POINTS = 300

class ASA(object):
    """
        Accessible surface area calculations
        """
    
    def __init__(self, pdb, unbound_a=None, unbound_b=None):
        """
        Calculates asa for pdb. you can specify pdbA and pdbB as unbound version
        otherwise it just extract it from the pdb
        """
        self.pdb = pdb
        if unbound_a is not None:
            self.pdbA = unbound_a[0]
            self.pdbB = unbound_b[0]
            self.chainA = unbound_a[1]
            self.chainB = unbound_b[1]
        else:
            self.pdbA = None
            self.pdbB = None
            self.chainA = None
            self.chainB = None
        self.interASA = dict()  # calculated from complex
        self.intraASA = dict()  # calculated from complex
        self.diffASAperAtom = dict()  # calculated from complex
        self.interPerAtom = dict()  # calculated from complex
        
        # not used now
        self.unboundASA = dict()  # unbound ASA per chain
        self.intraASAUnB = dict()  # calculated from unbound
        self.diffASAperAtomUnbound = dict()
        # results
        self.interfaceASA = 0  # ASA in interface bound against unbound
        self.combinedASA = 0  # TOTAL diff ASA (intra-inter from complex)
        self.boundVsUnboundASA = 0  # NOT USED
        self.unboundIntraTotal = 0  # NOT USED
        self.unboundMisses = 0
        # for outputing
        self.perAtom = self.pdb.getFile('per-atom.pdb')
        self.q_points = 0
        self.interfaceASAComp = {'A': 0, 'B': 0}
    
    def execute(self):
        """Calculate asa only for interface chains"""
        for chain in ''.join(self.pdb.interfaceParts):
            logging.info('calculating ASA for chain %s', chain)
            self.interASA[chain], self.intraASA[chain] = self.mesh(chain)
    
        intra_total = sum([sum([self.intraASA[chain] for chain in part]) for part in self.pdb.interfaceParts])
        inter_total = sum([sum([self.interASA[chain] for chain in part]) for part in self.pdb.interfaceParts])
        self.combinedASA = intra_total - inter_total
        
        if self.pdbA is not None and self.pdbB is not None:
            self.calc_interface()

#        self.pymol_script()
        self.perAtom.close()

    @staticmethod
    def calc_unbound_protein(pdb):
        """Calculates ASA for unbound protein
        :param pdb: protein
        :return: dictionary with ASA for each atom
        """
        asa_per_atom = dict()
        for chain in pdb.interfaceParts[0]:
            for atom, neighbors in ASA.nearAtomsOnChain(chain, pdb.atoms, pdb.ktree):
                intra_neighbors = [a for a in neighbors if a.pseudoChain == atom.pseudoChain]
                intra_area = ASA.calcASAforAtom(atom, intra_neighbors)
                asa_per_atom[atom] = intra_area

        return asa_per_atom

    def calc_interface(self):
        """Calculates ASA for interface"""
        diff_asa_a, misses_a = self.calcASAUnboundMissingReduced('A', self.pdbA, self.chainA)
        diff_asa_b, misses_b = self.calcASAUnboundMissingReduced('B', self.pdbB, self.chainB)
        
        self.interfaceASAComp['A'] = (diff_asa_a, misses_a)
        self.interfaceASAComp['B'] = (diff_asa_b, misses_b)
        self.interfaceASA = diff_asa_a + diff_asa_b
        self.unboundMisses = misses_a + misses_b
        logging.info('INTERFACE ASA: %s', str(self.interfaceASA))
    
    def calcASAUnboundMissingReduced(self, pseudo_chain, unbound_protein, unbound_chain):
        """calc interface ASA for bound-unbound
        interface: ASA=0 in inter
        :param pseudo_chain:
            :param unbound_protein:
            :param unbound_chain:
            :return:
            """
        ktree = unbound_protein.ktree
        inter_per_atom = self.interPerAtom
        
        def unbound_neighbors(atom):
            return [atom2 for atom2 in ktree.findByDistance(query_point=atom.coord, distance=(radio_atom(atom.atomType) + R_WATER * 2 + 1.8) ** 2) if atom2 != atom
                    and atom2.chain in unbound_chain]
        
        mapping = mapUnbound(self.pdb, unbound_protein, self.pdb.interfaceParts[0 if pseudo_chain == 'A' else 1],
                             unbound_chain)
        with self.pdb.getFile("missingASAInterface{0}.txt".format(unbound_protein.name)) as missInterface:
            diff_asa, misses = 0, 0
            for a_bound, asa in self.diffASAperAtom.items():
                if asa > 0 and a_bound.pseudoChain == pseudo_chain:
                    if a_bound not in mapping:
                        print('{0:^4}{1:^6}{2:^4}'.format(a_bound.chain, a_bound.resId, a_bound.symbol),
                              file=missInterface)
                        logging.warning('interface atom is not mapped %s', a_bound)
                        misses += 1
                        continue
                    a_unbound = mapping[a_bound]
                    neighbors = unbound_neighbors(a_unbound)
                    intra_area = ASA.calcASAforAtom(a_unbound, neighbors)
                    inter_area = inter_per_atom[a_bound]
                    diff_asa += (intra_area - inter_area)
                                                               
        logging.info('INTERFACE ASA: %s', str(self.interfaceASA))
        return diff_asa, misses
    
    def completeASA(self, missing_atoms):
        """Completes ASA for pdbUnbound based on missingAtoms
        How it works:
        On a specific residue, calculates the ASA from the first missing atom to
        the last atom in this residue
        The calculation is done on a residue alone (not dependent on neighbors residues)
        :param missing_atoms: list of atoms, from complex on the same chain
        :return: ASA to add the unbound
        """
        # first group missing atoms by res
        grouped = dict()
        for a in missing_atoms:
            if a.resId not in grouped:
                grouped[a.resId] = []
            grouped[a.resId].append(a)
        
        total_asa = 0
        pdb = self.pdb
        
        for resId, missAtoms in grouped.items():
            missAtoms.sort(key=lambda x: x.index)
            # first atom that appear iin both complex and unbound on the same residue
            first_non_missing = pdb[missAtoms[0].index - 1]
            
            tail = [a for a in pdb.atoms if
                    a.resId == resId and a.index > first_non_missing and a.chain == first_non_missing.chain]
            neighbors = ASA.nearAtomsOnChain(first_non_missing.chain, pdb.atoms, pdb.kdtree)
            neighbors = [a for a in neighbors if a.resId == resId and a.index < first_non_missing.index]
            # remove asa for the last atom in the residue
            total_asa -= ASA.calcASAforAtom(first_non_missing, neighbors)
            # add to asa area of the remaining tail
            for atom in tail:
                neighbors = [a for a in ASA.nearAtomsOnChain(atom.chain, pdb.atoms, pdb.kdtree) if a.resId == resId]
                total_asa += ASA.calcASAforAtom(atom, neighbors)
        return total_asa
    
    def pymol_script(self):
        """Generates pymol script for selection of the interface
        """
        for chain in ''.join(self.pdb.interfaceParts):
            logging.info('-------------')
            logging.info('ASA for chain %s', chain)
            logging.info('--inter: %s', self.interASA[chain])
            logging.info('--intra: %s', self.intraASA[chain])
            logging.info('diff ASA: %s', self.intraASA[chain] - self.interASA[chain])
        
        # 2 definitions for interface: based on distance and based on ASA
        interface_distance, interface_asa = dict(), dict()
        interfaces = self.pdb.getInterface()
        for atom in interfaces:
            if atom.chain not in interface_distance:
                interface_distance[atom.chain] = set()
            interface_distance[atom.chain].add(str(atom.resId))
        
        for atom in [a for a in self.pdb.atoms if (a in self.diffASAperAtom) and self.diffASAperAtom[a] > 0]:
            if atom.chain not in interface_asa:
                interface_asa[atom.chain] = set()
            interface_asa[atom.chain].add(str(atom.resId))
        
        sel_dist, sel_asa = [], []
        for chain, reses in interface_distance.items():
            sel_dist.append(' res ' + '+'.join(reses) + ' and chain ' + chain)
        for chain, reses in interface_asa.items():
            sel_asa.append(' res ' + '+'.join(reses) + ' and chain ' + chain)
        
        with self.pdb.getFile('.interface.pml') as interface_output:
            print('\n'.join(['sel ' + ' or '.join(sel_dist),
                             'show lines,sele',
                             'set_name sele,Distance-int',
                             'sel ' + ' or '.join(sel_asa),
                             'show lines,sele',
                             'set_name sele,ASA-int']), file=interface_output)

    def mesh(self, chain):
        """Create mesh of points around atoms in a chain
        :param chain:
        :return:
        """
        inter_asa, intra_asa = 0, 0
        interstring_chains = ''.join(self.pdb.interfaceParts)
        for atom, neighbors in ASA.nearAtomsOnChain(chain, self.pdb.atoms, self.pdb.ktree):
            # pseudo chain or chain?
            intra_neighbors = [a for a in neighbors if a.pseudoChain == atom.pseudoChain]
            inter_optimized = [a for a in neighbors if a.chain in interstring_chains]
            intraArea = ASA.calcASAforAtom(atom, intra_neighbors, callback=self.print_point)
            interArea = ASA.calcASAforAtom(atom, inter_optimized)
                    
            self.interPerAtom[atom] = interArea
            self.diffASAperAtom[atom] = intraArea - interArea
            inter_asa += interArea
            intra_asa += intraArea
        return inter_asa, intra_asa

    def getASA(self, atom):
        """Get ASA of an atom
        :param atom: atom to get ASA for
        :return: ASA of the atom
        """
        return self.diffASAperAtom[atom]

    def isBuried(self, atom):
        """Check whether atom is buried (inaccessible)
        :param atom: atom to check is buried
        :return: True is inaccessible/buried otherwise false
        """
        return self.interPerAtom[atom] == 0

    @staticmethod
    def calcASAforAtom(atom, neigbors, callback=None):
        r = radio_atom(atom.atomType) + R_WATER
        spiral_points = list(spiral(r, atom))
        
        for atom2 in neigbors:
            r2 = radio_atom(atom2.atomType)
            threshold = (r2 + R_WATER) ** 2
            for x, y, z in spiral_points[:]:
                dist = atom2.distanceFromXYZ((x, y, z))
                if dist < threshold:  # water is added to two sides of the equation?
                    spiral_points.remove((x, y, z))
        if callback:
            callback(spiral_points[:])
        
        area = area_calc(r, len(spiral_points), SPIRAL_POINTS)
        return area
    
    @staticmethod
    def nearAtomsOnChain(chain, atoms, ktree):
        extraRad = R_WATER * 2 + 1.8
        for atom in [atom for atom in atoms if atom.chain == chain]:
            neigbors = []
            for atom2 in ktree.findByDistance(query_point=atom.coord, distance=(radio_atom(atom.atomType) + extraRad) ** 2):
                if atom2 != atom:
                    neigbors.append(atom2)
            yield atom, neigbors

    def print_point(self, points):
        """Prints point to file
        :param points: points to print
        """
        for x, y, z in points:
            self.q_points += 1
            point_to_print = 'HETATM {0:^4}  O   HOH I {1:^4} {2:^11.3f} {3:^4.3f} {4:^8.3f}  0.83 56.58           O '.format(self.q_points, self.q_points, x, y, z)
            print(point_to_print, file=self.perAtom)


def radio_atom(atom):
    """Get the radius of an atom
        :param atom: atom to get radius for
        :return: radius of the atom
        """
    if atom == 'X':
        logging.error('undetermined atom while calculating ASA. Assuming C')
    return KNOWN_RADIUS[atom]


def area_calc(radius, point_in, total_points):
    """Calculates the partial area of ball
        :param radius: radius of ball
        :param point_in: points of the total points to include
        :param total_points: number of sampled points
        :return: area
        """
    return (4 * pi * radius ** 2) * point_in / total_points


def spiral(r, atom, n_points=SPIRAL_POINTS):
    """Generator for sampled points in ball
        :param r: radius of the ball
        :param atom: atom
        :param n_points: number of points
        """
    a = (4 * pi * r ** 2) / n_points
    d = sqrt(a)
    mPol = int(round((pi/d), 0))
    dPol = pi/mPol
    dAz = a/dPol
    for m in range(mPol):
        Pol = pi * (m + 0.5)/mPol
        mAz = int(round((2 * pi * sin(Pol)/dAz), 0))
        for n in range(mAz):
            Az = 2 * pi * n/mAz
            xCoord = (r * sin(Pol) * cos(Az)) + atom.x
            yCoord = (r * sin(Pol) * sin(Az)) + atom.y
            zCoord = (r * cos(Pol)) + atom.z
            yield (xCoord, yCoord, zCoord)


def main():
    """Script for calculating ASA"""
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('complex', help='PDB file of complex')
    parser.add_argument('complex_chains', nargs='+', help='Chain in the complex')
    parser.add_argument('ubnound_a', help='path to PDB of unbound component')
    parser.add_argument('ubnound_b', help='path to PDB of unbound component')
    parser.add_argument('chain_a', help='Chains in the first component')
    parser.add_argument('chain_b', help='Chains in the second component')
    args = parser.parse_args()
    
    print('ASA script')
    pdb = PDBReader.readFile(args.complex, interface_parts=args.complex_chains)
    pdbA = PDBReader.readFile(args.ubnound_a)
    pdbB = PDBReader.readFile(args.ubnound_b)
    asa = ASA(pdb, unbound_a=(pdbA, args.chain_a), unbound_b=(pdbB, args.chain_b))
    asa.execute()
    
    print(asa.combinedASA)
    print(asa.boundVsUnboundASA)


if __name__ == "__main__":
    logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)
    main()
