"""Hydrogen bonds module"""
from __future__ import print_function
import logging
import math
from math import sqrt

from .donorAcceptor import DonAcceptor
from .pdbReader import angle

MIN_LENGTH_HBOND = 1.97  # The typical length of a hydrogen bond in water is 197 pm
MAX_HBOND = 3.6 ** 2  # distance**2
HBOND_ENERGY = -0.01
WATER_LENGTH = 6
DEBUG = False


class hbonds(object):
    """Calculation and find of hydrogen bonds and hbonds with water on interface in complex.
    both intra-molecular  and inter-molecular are found
    """

    def __init__(self, pdb):
        self.pdb = pdb
        self.waterHbonds = 0  # number of hbonds with water on interface
        self.waterInterface = 0  # number of water moleculs on interface
        self.hbondsList = None
        self.HDPlusDefinition = True
        self.waterHbondsList = None

    def execute(self):
        interfaces = self.pdb.getInterface()
        # calculate water hbond
        self.h2Obonds(interfaces)
        # calculate hbonds
        self.hbondsList = self.hbonds(interfaces)

    def h2Obonds(self, interfaces):
        hbonds, minDistanceA, minDistanceB, interfaceB = self.buildWater(interfaces)
        self.H2Oprintout(hbonds, minDistanceA, minDistanceB, interfaceB)
        self.waterHbonds = len(hbonds)
        self.waterInterface = len(interfaceB)

    def buildWater(self, interfacesAAA):
        # H2O porting to python
        minDistanceA = dict()  # most closer atom from A to water
        minDistanceB = dict()  # same for B
        # water: (atomA, distance)
        for water in self.pdb.waters:
            minDistanceA[water] = (None, 7)
            minDistanceB[water] = (None, 7)

        probablInterface = set()
        lessThan35 = dict()  # Hcercano
        interfaceB = set()
        # foreach interface atom
        # for atom in [a for a in self.atoms if a.isInterface()]:
        for atom in [a for a in interfacesAAA if a.atomType != 'H']:
            for water in self.pdb.waters:
                sqDist = sqrt(water.distance(atom))
                if sqDist == 0:
                    raise Exception('WARNING: distance between water and atom is zero?')
                if sqDist < WATER_LENGTH:
                    if atom.pseudoChain == 'A':  # todo change the A as it could be B or E
                        probablInterface.add(water)  # (9)
                        if sqDist < minDistanceA[water][1]:
                            minDistanceA[water] = atom, sqDist
                    elif atom.pseudoChain == 'B' and (water in probablInterface) and min(minDistanceA[water][1],
                                                                                         sqDist) < 3.5:
                        if sqDist < minDistanceB[water][1]:
                            minDistanceB[water] = atom, sqDist
                        interfaceB.add(water)  # 11 in matrix is it same as 9
                    if sqDist < 3.5:
                        if water not in lessThan35:
                            lessThan35[water] = set()
                        lessThan35[water].add(atom)
        # atom water tuple list
        hbonds = []
        hbondMethod = self.waterHbondHDPlus if self.HDPlusDefinition else self.newWaterHbond
        for water in interfaceB:
            for atom in lessThan35[water]:
                dist = sqrt(water.distance(atom))
                donor, acceptor = hbondMethod(atom, water, dist)
                if donor or acceptor:
                    hbonds.append((atom, water))
        self.waterHbondsList = hbonds
        return hbonds, minDistanceA, minDistanceB, interfaceB

    def waterHbondHDPlus(self, atom, water, dist):
        """
            This method checks for hbond according to definition of hbplus
            see definition in: http://www.csb.yale.edu/userguides/datamanip/hbplus/hbplus_descrip.html
            this method doesn't give enegy just a boolean
        """
        MAX_DA = 3.9
        MAX_HA = 2.5
        MIN_DHA = 0.5 * math.pi
        MIN_D_A_AA = 0.5 * math.pi
        MIN_H_A_AA = 0.5 * math.pi
        if dist > MAX_DA:
            return (False, False)

        if atom.residue == 'HIS' and atom.symbol in ['ND1', 'NE2']:
            if dist > math.sqrt(2):
                return (True, True)

        for donor, donoarOrbital, hAtom, testD in DonAcceptor.donors(atom, self.pdb.getNextAtoms):
            DH = sqrt(donor.distance(hAtom))
            HA = sqrt(water.distance(hAtom))
            theta = angle(DH, dist, HA)
            # phi = 120/180*math.pi
            if HA < MAX_HA and (theta > MIN_DHA):  # and phi>MIN_H_A_AA
                return (True, True)

        for acceptor, acceptorOrbital, preAtom, testA in DonAcceptor.acceptors(atom, self.pdb.getNextAtoms):
            # theta = math.pi #180
            phi = angle(sqrt(preAtom.distance(acceptor)), sqrt(preAtom.distance(water)), dist)
            preAtomAcceptor = sqrt(preAtom.distance(acceptor))
            preAtomDonor = sqrt(preAtom.distance(water))
            d_a_aa = angle(dist, preAtomDonor, preAtomAcceptor)
            # 3.5 = 2.5 (MAX_HA)+1A from water
            if dist < 3.5 and phi > MIN_H_A_AA and d_a_aa > MIN_D_A_AA:  # and theta>MIN_DHA and HA<MAX_HA
                return (True, True)

        return (False, False)

    def newWaterHbond(self, atom, water, dist):
        haveHbond = False
        eHB = 0
        if dist < MIN_LENGTH_HBOND:
            logging.info('WARNING no hbond between', atom, water, '. distance is too small: {0:.3f}'.format(dist))
            return False, 0
        for donor, donoarOrbital, hAtom, testD in DonAcceptor.donors(atom, self.pdb.getNextAtoms):
            acceptorOrbital = 'sp3'  # water is sp3
            DH = sqrt(atom.distance(hAtom))
            HA = sqrt(water.distance(hAtom))
            theta = angle(DH, dist, HA)
            phi = 120 / 180 * math.pi
            f = self.fAngleDep(donoarOrbital, acceptorOrbital, theta, phi)
            energy = self.leenardJones(dist) * f
            if energy < HBOND_ENERGY:
                if haveHbond:
                    eHB = min(eHB, energy)
                else:
                    eHB = energy
                haveHbond = True

        for acceptor, acceptorOrbital, preAtom, testA in DonAcceptor.acceptors(atom, self.pdb.getNextAtoms):
            donoarOrbital = 'sp3'  # water is sp3
            theta = math.pi  # 180
            phi = angle(sqrt(preAtom.distance(acceptor)), sqrt(preAtom.distance(water)), dist)
            f = self.fAngleDep(donoarOrbital, acceptorOrbital, theta, phi)
            energy = self.leenardJones(dist) * f
            if energy < HBOND_ENERGY:
                if haveHbond:
                    eHB = min(eHB, energy)
                else:
                    eHB = energy
                haveHbond = True
        return haveHbond, haveHbond

    def hbonds(self, interfaces):
        """
        used by hbonds not by hbonds water
        """

        logging.debug('start hbonds for interface atoms')
        donorsAcceptors = set()
        for atom in [a for a in interfaces]:
            # dont ignore by chain - both inter and intra
            for atom2 in [a for a in self.pdb.ktree.findByDistance(query_point=atom.coord, distance=MAX_HBOND) if
                          a != atom]:
                # when atom/atom2 can be donor/acceptor or acceptor/donor
                hbondMethod = self.checkHbondHBPlus if self.HDPlusDefinition else self.checkHbond
                hbond, eHb = hbondMethod(atom, atom2)
                if hbond:
                    donorsAcceptors.add((atom, atom2, eHb))
                hbond, eHb = hbondMethod(atom2, atom)
                if hbond:
                    donorsAcceptors.add((atom2, atom, eHb))

        self.hbondsOutput(donorsAcceptors)
        logging.info('Found %s hbonds', len(donorsAcceptors))
        return donorsAcceptors

    def fAngleDep(self, dOrbital, aOrbital, theta, phi):
        PhiFixRad = 109.5 / 180 * math.pi
        Rad90 = 0.5 * math.pi

        dOrbital = dOrbital.strip()
        aOrbital = aOrbital.strip()
        if dOrbital == 'sp3' and aOrbital == 'sp3':
            if theta > Rad90 and phi - PhiFixRad < Rad90:
                f = (math.cos(theta) ** 2) * (math.cos(phi - PhiFixRad) ** 2)
            else:
                logging.info('sp3-sp3 criteria failed')
                f = 0
        elif dOrbital == 'sp3' and aOrbital == 'sp2':
            if phi > Rad90:
                f = (math.cos(theta) ** 2) * (math.cos(phi) ** 2)
            else:
                logging.info('sp3-sp2 criteria failed')
                f = 0
        elif dOrbital == 'sp2' and aOrbital == 'sp3':
            if theta > Rad90:  # our new condition not in article
                f = (math.cos(theta) ** 4)
            else:
                logging.info('sp2-sp3 criteria failed')
                f = 0
        elif dOrbital == 'sp2' and aOrbital == 'sp2':
            if theta > Rad90:  # our new condition not in article
                f = (math.cos(theta) ** 2) * (math.cos(phi) ** 2)
            else:
                logging.info('sp2-sp2 criteria failed')
                f = 0
        else:
            raise Exception('error dOrbital:(' + dOrbital + ')  aOrbital:(', aOrbital, ')')
        return f

    def leenardJones(self, r):
        # 12-10 LJ
        D0 = 8
        R0 = 2.8
        return D0 * (5 * (R0 / r) ** 12 - 6 * (R0 / r) ** 10)

    def checkHbondHBPlus(self, pDonor, pAcceptor):
        MAX_DA = 3.9
        MAX_HA = 2.5
        MIN_DHA = 0.5 * math.pi
        MIN_D_A_AA = 0.5 * math.pi
        MIN_H_A_AA = 0.5 * math.pi

        """
            This method checks for hbond according to definition of hbplus
            see definition in: http://www.csb.yale.edu/userguides/datamanip/hbplus/hbplus_descrip.html
            this method doesn't give enegy just a boolean
        """
        if pDonor.residue == 'HIS' and pDonor.symbol in ['ND1', 'NE2']:
            for acceptor, acceptorOrbital, preAtom, testA in DonAcceptor.acceptors(pAcceptor, self.pdb.getNextAtoms):
                dist = sqrt(pAcceptor.distance(pDonor))
                preAtomAcceptor = sqrt(preAtom.distance(acceptor))
                preAtomDonor = sqrt(preAtom.distance(pDonor))
                d_a_aa = angle(dist, preAtomDonor, preAtomAcceptor)
                if dist > math.sqrt(2) and dist < 3.9 and d_a_aa > MIN_D_A_AA:
                    return (True, 0)

        for donor, donoarOrbital, hAtom, testD in DonAcceptor.donors(pDonor, self.pdb.getNextAtoms):
            for acceptor, acceptorOrbital, preAtom, testA in DonAcceptor.acceptors(pAcceptor, self.pdb.getNextAtoms):
                DH = sqrt(donor.distance(hAtom))
                HA = sqrt(pAcceptor.distance(hAtom))
                dist = sqrt(pAcceptor.distance(donor))
                if dist < MAX_DA and HA < MAX_HA:
                    theta = angle(DH, dist, HA, radians=True)
                    preAtomAcceptor = sqrt(preAtom.distance(acceptor))
                    preAtomDonor = sqrt(preAtom.distance(donor))
                    phi = angle(preAtomAcceptor, sqrt(preAtom.distance(hAtom)), sqrt(acceptor.distance(hAtom)),
                                radians=True)
                    d_a_aa = angle(dist, preAtomDonor, preAtomAcceptor)
                    # for dbeug
                    # return (True, 0)
                    if (theta > MIN_DHA and phi > MIN_H_A_AA and d_a_aa > MIN_D_A_AA):
                        return (True, 0)
        return (False, 0)

    def checkHbond(self, pDonor, pAcceptor):
        """
            This is extension to waterHbond
            atom is the donor and atom2 is the acceptor
            This method is based on leenard jones and gives energy pontencial for each possible bond
        """
        global DEBUG
        check = False
        haveHbond = False
        eHB = 0
        # if the atom is the donor and water is the acceptor
        for donor, donoarOrbital, hAtom, testD in DonAcceptor.donors(pDonor, self.pdb.getNextAtoms):
            for acceptor, acceptorOrbital, preAtom, testA in DonAcceptor.acceptors(pAcceptor, self.pdb.getNextAtoms):
                DH = sqrt(donor.distance(hAtom))
                HA = sqrt(pAcceptor.distance(hAtom))
                dist = sqrt(pAcceptor.distance(donor))
                if dist < MIN_LENGTH_HBOND:
                    logging.info('WARNING no hbond between', donor, pAcceptor,
                                 '. distance is too small: {0:.3f}'.format(dist))
                    return False, 0
                theta = angle(DH, dist, HA, radians=True)
                phi = angle(sqrt(preAtom.distance(acceptor)), sqrt(preAtom.distance(hAtom)),
                            sqrt(acceptor.distance(hAtom)), radians=True)
                logging.info(
                    'Orbitals Donor: %s Acceptor: %s . Theta %s Phi %s' % (donoarOrbital, acceptorOrbital, theta, phi))
                f = self.fAngleDep(donoarOrbital, acceptorOrbital, theta, phi)
                energy = self.leenardJones(dist) * f
                if DEBUG:
                    # if donor.resId in [63,23]:
                    print('%s%s%s-%s%s%s E=%s Dist=%s Phi=%s Theta=%s ' % (
                    donor.chain, donor.resId, donor.symbol, acceptor.chain, acceptor.resId, acceptor.symbol, energy,
                    dist, phi, theta))
                if energy < HBOND_ENERGY:
                    if haveHbond:
                        eHB = min(eHB, energy)
                    else:
                        eHB = energy

                    haveHbond = True
        return (haveHbond, eHB)

    def hbondsOutput(self, donorsAcceptors):

        with self.pdb.getFile('.hbonds.txt') as hBondsoutput:
            print("			 " + self.pdb.name + "\n", file=hBondsoutput)
            chainDonors = dict()
            chainAcceptors = dict()
            for donor, acceptor, eHb in donorsAcceptors:
                if donor.chain not in chainDonors:
                    chainDonors[donor.chain] = set()
                if acceptor.chain not in chainAcceptors:
                    chainAcceptors[acceptor.chain] = set()
                chainDonors[donor.chain].add(str(donor.resId))
                chainAcceptors[acceptor.chain].add(str(acceptor.resId))
                print('{0:<6}{1:<3}{2:<3} (Donor) <--> {3:<6} {4:<3}  {5:<3} (Acceptor)   Energy = {6}'.format(
                    donor.symbol, donor.chain, donor.resId, acceptor.symbol, acceptor.chain, acceptor.resId, eHb),
                    file=hBondsoutput)

        # TODO: if atom is acceptor to many
        with self.pdb.getFile('.hbonds.pml') as pmlHbonds:
            donorSelection, accSelection = [], []
            for chain, donorsList in chainDonors.items():
                donorSelection.append(' res ' + '+'.join(donorsList) + ' and chain ' + chain)
            for chain, acceptorList in chainAcceptors.items():
                accSelection.append(' res ' + '+'.join(acceptorList) + ' and chain ' + chain)

            pml_script = '\n'.join([
                'sel ' + ' or '.join(donorSelection),
                'show lines,sele',
                'set_name sele,Donor',
                'sel ' + ' or '.join(accSelection),
                'show lines,sele',
                'set_name sele,Acceptor'
            ])
            print(pml_script, file=pmlHbonds)

    def H2Oprintout(self, hbonds, minDistanceA, minDistanceB, interfaceB):
        # h2o interface
        with self.pdb.getFile('.interfaces.txt') as interfaceFile:
            print("			 " + self.pdb.name + "\n", file=interfaceFile)

            for water in self.pdb.waters:
                if water in interfaceB:
                    atomA, neighborDis = minDistanceA[water]
                    atomB, neighborDis = minDistanceB[water]
                    strLine = "{0:^3} {1}   -   {2:^4} {3:^3} {4}  {5:^4}  & {6:^4} {7:^3} {8}   {9:^4}".format(
                        water.resId, water.chain,
                        atomA.symbol, atomA.residue, atomA.chain, atomA.resId,
                        atomB.symbol, atomB.residue, atomB.chain, atomB.resId)
                    print(strLine, file=interfaceFile)

        # hbonds file
        with self.pdb.getFile('.h2oHbonds.txt') as hbondsFile:
            print("			 " + self.pdb.name, file=hbondsFile)
            print("	      (" + str(len(hbonds)) + " Hbonds) \n", file=hbondsFile)

            for atom, water in hbonds:
                print("HOH  {0:^2} {1:^2}   -   {2:^3} {3:^4} {4:^2} {5:^4}".format(water.chain, water.resId,
                                                                                    atom.symbol, atom.residue,
                                                                                    atom.chain, atom.resId),
                      file=hbondsFile)
            hbondsFile.close()
            print('found  %i hbonds with water on interface' % len(hbonds))
