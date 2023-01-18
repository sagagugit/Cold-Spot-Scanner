"""
Calculates VDW:
* according to drieding
* with cutoff of 9A
* between atoms in interface of two parts
* interface definition according to both asa and distance 4
"""
import math
import sys

sys.path.append('../')
from ..pdbReader import PDBReader
from ..ASA import radio_atom, spiral, KNOWN_RADIUS
from ..kdtree import KDTree
from .. import DBConfig

MAX_WW_RAD = max(KNOWN_RADIUS.values())


def getHbondsH(pdb):
    conn = DBConfig.get_connection()
    cursor = conn.cursor()

    cursor.execute("""select DonorChain,DonorResId,Hydrogen from
                        Ndrieding
                        inner join donors2
                        on DonorSymbol=donors2.Symbol
                        where PDB='%s'""" % pdb.name)
    interfaceAtoms = []
    for chain, resId, symbol in cursor.fetchall():
        interfaceAtoms += [a for a in pdb.atoms if a.chain == chain and a.resId == resId and a.symbol == symbol]
    return interfaceAtoms


def get_interface_atoms(pdb):
    conn = DBConfig.get_connection()
    cursor = conn.cursor()

    #get interface table
    cursor.execute("""select Chain,ResId,Symbol from
                        NinterfaceAtoms
                        where PDB='%s'""" % ( pdb.name ))

    interfaceAtoms = []
    for chain, resId, symbol in cursor.fetchall():
        interfaceAtoms += [a for a in pdb.atoms if a.chain == chain and a.resId == resId and a.symbol[0:3] == symbol]
    return interfaceAtoms


def getNamesAndChainFromDB():
    conn = DBConfig.get_connection()
    cursor = conn.cursor()

    #create interface table
    cursor.execute("""select Complex, niceLog from origTable2""")
    res = []
    for c, lg in cursor.fetchall():
        res.append((c[0:4], c[5:].split(':'), lg))
    return res


def calcFunc(atomA, atomBs, func, estimate):
    atomRadA = radio_atom(atomA.atomType)
    if len(atomBs) == 1:#exact
        atomB = atomBs[0]
        dist = math.sqrt(atomA.distance(atomB))
        atomRadB = radio_atom(atomB.atomType)
        collosionRad = atomRadB + atomRadA - dist
        collosionRad2 = (collosionRad / 2) ** 2

        baseIntA = math.sqrt(atomRadA ** 2 - collosionRad2)
        baseIntB = math.sqrt(atomRadB ** 2 - collosionRad2)
        #integral of ball surface according to radius
        volume1 = func(baseIntA)
        volume2 = func(baseIntB)
        return volume1 + volume2
    else:#estimate
        ballPoints = list(spiral(radio_atom(atomA.atomType), atomA))
        clashs = []

        for b in atomBs:
            for p in ballPoints:
                if p not in clashs and b.distanceFromXYZ(p) < radio_atom(b.atomType) ** 2:
                    clashs.append(p)
        volume = estimate(atomRadA)
        volume = volume * float(len(clashs) / len(ballPoints))
        return volume


def calcClash(atomA, atomBs, surface):
    if len(atomBs) == 0:
        return 0.0
    func = (lambda r: math.pi * r ** 2) if surface else (lambda r: math.pi * (r ** 3) / 3.0)
    #the estimator for volume assume the volume of 1/2 radius
    estimateFunc = (lambda r: 4.0 * math.pi * r ** 2) if surface else (lambda r: math.pi * ((r / 2.0) ** 3) * 4.0 / 3.0)
    return calcFunc(atomA, atomBs, func, estimateFunc)


def calcLJ(x, y, dist, hHbonds):
    #see p3 from DRIEDING: A generic force field for molecule simulations
    R0dic = {
        'H': 3.195,
        'C': 3.8983,
        'N': 3.6621,
        'O': 3.4046,
        'S': 4.03
    }
    D0Dict = {
        'H': 0.0152,
        'C': 0.0951,
        'N': 0.0774,
        'O': 0.0957,
        'S': 0.3440
    }

    D0X = D0Dict[x.atomType]
    D0Y = D0Dict[y.atomType]
    if x in hHbonds:
        D0X = 0.0001
        #return 0.0
    if y in hHbonds:
        D0Y = 0.0001
        #return 0.0
    D0 = math.sqrt((D0X * D0Y))
    R0 = math.sqrt((R0dic[x.atomType] * R0dic[y.atomType]))
    #use the same as AMBER and DRIEDING
    #use aritmethic mean instead of geometric mean can be used instead
    #"more consistent with chemicel practice"
    p = (dist / R0)
    return D0 * (p ** -12 - 2 * p ** -6)


def calcX6(x, y, dist):
    M0Dict = {
        'H': 12.382,
        'C': 14.034,
        'N': 13.843,
        'O': 13.483,
        'S': 12.0
    }
    D0Dict = {
        'H': 0.0152,
        'C': 0.0951,
        'N': 0.0774,
        'O': 0.0957,
        'S': 0.3440
    }
    R0dic = {
        'H': 3.195,
        'C': 3.8983,
        'N': 3.6621,
        'O': 3.4046,
        'S': 4.03
    }
    D0 = math.sqrt((D0Dict[x.atomType] * D0Dict[y.atomType]))
    R0 = math.sqrt((R0dic[x.atomType] * R0dic[y.atomType]))
    p = (dist / R0)
    M0 = 13.772
    return D0 * ((6.0 / (M0 - 6)) * math.exp(M0 * (1 - p)) - (M0 / (M0 - 6)) * p ** -6)


def calcWDW(atomA, atomBs, lj=True, hHbonds=[]):
    sumE = 0.0
    for b in atomBs:
        distance = math.sqrt(b.distance(atomA))
        e = calcLJ(atomA, b, distance, hHbonds) if lj else calcX6(atomA, b, distance)
        sumE += e

    return sumE


def calcCompl(pdb_path, chains):
    pdb = PDBReader.readFile(pdb_path, chains)

    interface = get_interface_atoms(pdb)
    hHbonds = getHbondsH(pdb)
    components = []
    for part in pdb.interfaceParts:
        components.append([a for a in interface if a.chain in part])

    comp2 = [atom for atom in pdb.atoms if atom.chain in pdb.interfaceParts[1]]
    ktree = KDTree.construct_from_data(comp2)
    clashV, clashS = 0.0, 0.0
    sumVDW = 0.0
    sumVDWx = 0.0

    for atom in components[0]:
        atomRad = radio_atom(atom.atomType)
        nearAtoms = list(ktree.findByDistance(query_point=atom.coord, distance=9 ** 2))#(MAX_WW_RAD+atomRad+EPSILON)
        collusion = [a for a in nearAtoms if math.sqrt(atom.distance(a)) < radio_atom(a.atomType) + atomRad]
        contact = [a for a in nearAtoms if a not in collusion]
        surfaceRep = calcClash(atom, collusion, True)
        volumeRep = calcClash(atom, collusion, False)
        vdw = calcWDW(atom, contact, hHbonds=hHbonds)
        vdwx6 = calcWDW(atom, contact, lj=False)
        sumVDW += vdw
        sumVDWx += vdwx6
        clashS += surfaceRep
        clashV += volumeRep

    return sumVDW, sumVDWx, clashV, clashS
