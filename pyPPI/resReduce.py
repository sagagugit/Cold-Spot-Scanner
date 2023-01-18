import logging

from .alignment.NWAlignment import NWAlignment
from .atom import ResiduesCodes

PRINT_ALIGNED = False


def getResidues(pdbBound, pdbUnbound, chainBound, chainUnbound):
    resBound, resUnbound = [], []
    for chain in chainBound:
        resBound += [ResiduesCodes[a.residue if a.residue != 'GLX' else 'GLN'] for a in pdbBound.atoms if
                     a.chain == chain and a.symbol == 'CA']
    for chain in chainUnbound:
        resUnbound += [ResiduesCodes[a.residue if a.residue != 'GLX' else 'GLN'] for a in pdbUnbound.atoms if
                       a.chain == chain and a.symbol == 'CA']
    return resBound, resUnbound


def getAtoms(pdbBound, pdbUnbound, chainBound, chainUnbound):
    resBound, resUnbound = [], []
    for chain in chainBound:
        resBound += [a for a in pdbBound.atoms if a.chain == chain]
    for chain in chainUnbound:
        resUnbound += [a for a in pdbUnbound.atoms if a.chain == chain]
    return resBound, resUnbound


def alignResidues(bound, unbound):
    global PRINT_ALIGNED
    seq1 = ''.join(bound)
    seq2 = ''.join(unbound)
    align = NWAlignment(seq1, seq2)
    align.fillIn()
    nSeq1, nSeq2 = align.getTraceback()
    if PRINT_ALIGNED:
        printAligned(''.join(nSeq1), ''.join(nSeq2))

    resIdA, resIdB = -1, -1
    missesA, missesB = [], []
    for a, b in zip(nSeq1, nSeq2):
        if a != '-':
            resIdA += 1
        if b != '-':
            resIdB += 1

        if a == '-':
            missesB.append(resIdB)
        if b == '-':
            missesA.append(resIdA)

    return missesA, missesB


def printAligned(seqA, seqB):
    linesLength = 100
    for i in range(0, min(len(seqA), len(seqB)), linesLength):
        print('A  ' + seqA[i:i + linesLength])
        print('B  ' + seqB[i:i + linesLength])


def redcueAtoms(pdbBound, pdbUnbound, chainBound, chainUnbound):
    """
    Returns atoms of bound align to unbound
    When there is a miss within the range it is removed
    """
    resBound, resUnbound = getResidues(pdbBound, pdbUnbound, chainBound, chainUnbound)
    missesBound, missesUnbound = alignResidues(resBound, resUnbound)
    boundAtoms, unboundAtoms = getAtoms(pdbBound, pdbUnbound, chainBound, chainUnbound)

    mappingB = dict(((a.chain, a.resId), i) for i, a in enumerate([a for a in boundAtoms if a.symbol == 'CA']))
    boundReduce = [a for a in boundAtoms if
                   (a.chain, a.resId) in mappingB and mappingB[(a.chain, a.resId)] not in missesBound]
    for a in boundAtoms:
        if (a.chain, a.resId) not in mappingB:
            logging.warning('WARNING missing CA in ', a)
    mappingUn = dict(((a.chain, a.resId), i) for i, a in enumerate([a for a in unboundAtoms if a.symbol == 'CA']))
    unboundReduce = [a for a in unboundAtoms if
                     (a.chain, a.resId) in mappingUn and mappingUn[(a.chain, a.resId)] not in missesUnbound]
    for a in unboundAtoms:
        if (a.chain, a.resId) not in mappingUn:
            logging.warning('WARNING missing CA in ', a)
    return boundReduce, unboundReduce


def mapUnbound(pdbBound, pdbUnbound, chainBound, chainUnbound):
    """
        return a dictionary mapping between atom bound to atom unbound
    """
    boundReduce, unboundReduce = redcueAtoms(pdbBound, pdbUnbound, chainBound, chainUnbound)
    mapping = dict()
    unboundReduce.append(None)
    unboundIter = iter(unboundReduce)
    boundRes, unBoundRes = -1, -1
    residueUnbound = []
    y = next(unboundIter)
    for x in boundReduce:
        if boundRes != x.resId and y is not None:
            boundRes = x.resId
            residueUnbound = [y]
            unBoundRes = y.resId
            y = next(unboundIter)
            while y is not None and y.resId == unBoundRes:
                residueUnbound.append(y)
                y = next(unboundIter)

        for aUnbound in residueUnbound:
            if aUnbound.symbol == x.symbol:
                mapping[x] = aUnbound
    return mapping
