"""
Finds depth and periphriality for atoms
"""

import math

import numpy as np

from .. import DBConfig
from ..ASA import ASA
from ..kdtree import KDTree
from ..pdbReader import PDBReader


def get_interface_atoms(pdb, includingDistance=False):
    conn = DBConfig.get_connection()
    cursor = conn.cursor()

    if includingDistance:
        cursor.execute("""select Chain,ResId,Symbol from
                             NinterfaceAtoms
                             where PDB='%s'""" % pdb.name)
    else:
        cursor.execute("""select Chain,ResId,Symbol from
                             interfaceAtoms
                             where PDB='%s'""" % pdb.name)
    # find unbounds relating to
    interfaceAtoms = []
    for chain, resId, symbol in cursor.fetchall():
        interfaceAtoms += [a for a in pdb.atoms if a.chain == chain and a.resId == resId and a.symbol[0:3] == symbol]
    return interfaceAtoms


def assign_depth(interface):
    """
     Finds the distance of farest atom in interface for each atom
     """

    from scipy.spatial.distance import pdist, squareform
    interface_coords = np.array([a.coord for a in interface])
    depth = squareform(pdist(interface_coords)).max(1)
    maxDepth = depth.max()
    res = []

    for atom, dist in zip(interface, depth/maxDepth):
        res.append((atom, dist, dist / maxDepth))
    return res


def assign_peripheral(interface, surface_atoms):
    """
     Finds the nearest surface atom which is non interface to each atom
     """
    peripheral = []
    res = []
    aTree = KDTree.construct_from_data([atom for atom in surface_atoms if atom not in interface])
    max_peripheral = 0
    for atom in interface:
        near, distance2 = aTree.findNearest(query_point=atom.coord, num=1)
        max_peripheral = max(max_peripheral, distance2)
        peripheral.append((atom, distance2))
    max_peripheral = math.sqrt(max_peripheral)
    for atom, distance2 in peripheral:
        dist = math.sqrt(distance2)
        res.append((atom, dist, dist / max_peripheral))
    return res


def calc_peripheral_PDB(pdb_path, chains):
    pdb = PDBReader.readFile(pdb_path, chains)
    asa = ASA(pdb)
    asa.execute()
    interface = get_interface_atoms(pdb, True)
    components = list(interface)
    surfaceComponents = [atom for atom in asa.interPerAtom.keys() if not asa.isBuried(atom)]
    depth = assign_depth(components)
    peripheral = assign_peripheral(components, surfaceComponents)
    return depth, peripheral
