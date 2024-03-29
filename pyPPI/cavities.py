"""
Monte carlo estimation of cavity volume.
"""

import numpy as np
import os

from scipy.linalg import norm
from scipy.spatial import cKDTree

from Bio.PDB.PDBParser import PDBParser
import Bio.PDB.ResidueDepth as ResidueDepth
#from pyPPI.ResidueDepthcopy import get_surface, min_dist
from pyPPI.ResidueDepthcopy import get_surface

from .ASA import KNOWN_RADIUS, R_WATER
from .DBConfig import get_connection
from Bio.PDB import Select, PDBIO
from Bio.PDB.PDBParser import PDBParser

RESULTS_DIR = './results/'

class ChainSelect(Select):
    def __init__(self, chain):
        self.chain = chain

    def accept_chain(self, chain):
        if chain.get_id() in self.chain:
            return 1
        else:          
            return 0

def print_rand_points(pdb, points):
    """Prints rand point to file
        :param points: points to print
        """
    cavity_file =  pdb.name + '_cavities.pdb'
    with open(cavity_file, 'w') as rand_point_file:
        q_points = 0
        res = 1
        atoms_per_res = find_atoms_per_res(len(points))
        for x, y, z in points:
            q_points += 1
            point_to_print = 'HETATM{0:^5}  O   HOH I{1:^4}    {2:^7.3f} {3:^7.3f} {4:^7.3f}  1.00 43.38          O'.format(q_points, res, x, y, z)
            print(point_to_print, file=rand_point_file)
            if atoms_per_res == 0:
                res += 1
            else:
                if q_points % (10 * atoms_per_res) == 0:
                    res +=1
    
def find_atoms_per_res(number_of_atoms):
    if number_of_atoms < 10000 and number_of_atoms > 9999:
        return 1
    elif len(str(int(number_of_atoms))) <= 4:
        return 0
    else:
        return 1 + find_atoms_per_res(number_of_atoms / 10)

def calculateVolume(pdb, interface):
    
    rand_points, pointsToCheck, allVolume = createRandPoints(interface)
    print('Number of total points:', pointsToCheck)
    rand_points = inInterface(rand_points, pdb, interface)
    print('Points in Interface:', len(rand_points))
    rand_points = onlyBuried(rand_points, pdb)
    print('Points Buried:', len(rand_points))
    rand_points = notinVDWspheres(rand_points, pdb)
    print('Points out of van der Waals spheres:', len(rand_points))
    
    print_rand_points(pdb, rand_points)
    
    cavities_candidates = len(rand_points)
    cavitiesVolume = (cavities_candidates / pointsToCheck) * allVolume

    return cavitiesVolume
    
def createRandPoints(interface, top_reduction_factor=0.05, bottom_reduction_factor=0.01):
    ACCURACY = 1.0 ** 3
    ACCURACY_FACTOR = 10.0

    # Find bounding box with reduced top and bottom size
    boundingBox = []
    for i in range(0, 3):
        axisCoords = [a.coord[i] for a in interface]
        minAxis = min(axisCoords)
        maxAxis = max(axisCoords)

        # Calculate the reduced size for top and bottom
        top_reduction = (maxAxis - minAxis) * top_reduction_factor
        bottom_reduction = (maxAxis - minAxis) * bottom_reduction_factor

        # Apply reductions to top and bottom bounds
        boundingBox.append((minAxis + bottom_reduction, maxAxis - top_reduction))

    # Convert boundingBox to NumPy array
    boundingBox = np.array(boundingBox)

    distances = boundingBox[:, 1] - boundingBox[:, 0]
    allVolume = np.prod(distances)
    pointsToCheck = int((ACCURACY_FACTOR * (allVolume / ACCURACY)) * 2)
    np.random.seed(2023)
    rand_points = np.random.rand(pointsToCheck, 3)
    for i, min_max in enumerate(boundingBox):
        rand_points[:, i] *= min_max[1] - min_max[0]
        rand_points[:, i] += min_max[0]

    return rand_points, pointsToCheck, allVolume

def inInterface(rand_points, pdb, interface):
    
    components = [[a for a in interface if a.chain in part] for part in pdb.interfaceParts]
    partResIds = [set((a.chain, a.resId, a.iCode) for a in part) for part in components]
    interfaceRes = [[a for a in pdb.atoms if (a.chain, a.resId, a.iCode) in partResIds] for part, partResIds
                    in zip(pdb.interfaceParts, partResIds)]
    
    aTree = cKDTree(np.array([a.coord for a in interfaceRes[0]]))
    bTree = cKDTree(np.array([a.coord for a in interfaceRes[1]]))
    
    aDistances, neighbors_a = aTree.query(rand_points)
    bDistances, neighbors_b = bTree.query(rand_points)
    
    nearA = np.array([a for a in interfaceRes[0]])[neighbors_a]
    nearB = np.array([a for a in interfaceRes[1]])[neighbors_b]

    ab_Distance = norm(np.array([a.coord for a in nearA]) - np.array([a.coord for a in nearB]), axis=1)
    
    selector = (aDistances < ab_Distance) & (bDistances < ab_Distance)
   
    return rand_points[selector]


def onlyBuried(rand_points, pdb):
    
    parser = PDBParser(PERMISSIVE=1)
#    structure = parser.get_structure(pdb.name, pdb.file)
#    chains = ['A', 'B', 'C']
#    pdb_chain_file = '%s_C.pdb' % pdb.name                          
#    x = PDBIO()               
#    x.set_structure(structure)
#    x.save('{}'.format(pdb_chain_file), ChainSelect(chains))
    stru = parser.get_structure(pdb.name, '%s_C.pdb' % pdb.name )
    model = stru[0]
    surface = get_surface(model, probe_radius = 2.4)
    debug_surface(surface, pdb)
    
    pdbTree = pdb.ktree
    _, atomDistances = pdbTree.findNearest(query_point=rand_points, num=1)
    
    surfaceTree = cKDTree(surface)
    surfaceDistances, _ = surfaceTree.query(rand_points)
    
    selector = (atomDistances < surfaceDistances) & (surfaceDistances > 1.0)
    
    rand_points = rand_points[selector]
    return rand_points

def debug_surface(surface, pdb):
    debug_file = pdb.getFile('surface_points.pdb')
    atoms_per_res = find_atoms_per_res(len(surface))
    q_points = 0
    res = 1
    for x, y, z in surface:
        q_points += 1
        point_to_print = 'HETATM{0:^5}  O   HOH I{1:^4}    {2:^7.3f} {3:^7.3f} {4:^7.3f}  1.00 43.38          O'.format(q_points, res, x, y, z)
        print(point_to_print, file=debug_file)
        if atoms_per_res == 0:
            res += 1
        else:
            if q_points % (10 * atoms_per_res) == 0:
                res +=1
    debug_file.close()

def notinVDWspheres(rand_points, pdb):
    for element in KNOWN_RADIUS.keys():
        element_coords = [a.coord for a in pdb.atoms if a.atomType == element]
        if not any(element_coords):
            continue
        elementTree = cKDTree(element_coords)
        points_in_VDW_radius = elementTree.query_ball_point(rand_points, KNOWN_RADIUS[element] + R_WATER)
        selector = (points_in_VDW_radius.astype(bool) == False)
        rand_points = rand_points[selector]
        if not np.any(rand_points):
            print('No cavities in ' + pdb.name)
            break
    
    return rand_points
