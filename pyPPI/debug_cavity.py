import os
import numpy as np
from scipy.spatial import cKDTree
from scipy.linalg import norm
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.ResidueDepth import get_surface

def read_file(pdb_debug_file):
    rand_points = []
    with open(pdb_debug_file, 'r') as debug_file:
        for line in debug_file.readlines():
            x = float(line[30:38].rstrip()); y = float(line[38:46].rstrip()); z = float(line[46:54].rstrip())
            rand_points.append([x, y, z])
    return np.array(rand_points)

def debugging(pdb, interface):
    for pdb_debug_files in os.listdir('debug_rand_points'):
        if pdb_debug_files[:4] == pdb.name:
            pdb_debug_file = os.path.join('debug_rand_points', pdb_debug_files)
            break
    rand_points = read_file(pdb_debug_file)
    interfaceRes = findInterfaceRes(pdb, interface)
    SurfaceTree = get_surfaceTree(pdb)
    with open('debugging_file.txt', 'w') as debug_file:
        print('Debugging File for 3MTN mutant pdb', file=debug_file)
        surface_debugging(SurfaceTree, rand_points, debug_file)
        interface_debugging(interfaceRes, rand_points, debug_file)

def get_surfaceTree(pdb):
    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure(pdb.name, pdb.file)
    model = structure[0]
    surface = get_surface(model)
    SurfaceTree = cKDTree(surface)
    return SurfaceTree

def surface_debugging(SurfaceTree, rand_points, debug_file):
    surface_debug_points = [907, 2106]
    for point in surface_debug_points:
        print('\nPoints too clost to the surface', file=debug_file)
        print(('\nRand_point '+ str(point) + ': ' + rand_point[point - 1]), file=debug_file)
        depth, _ = SurfaceTree.query(rand_point[point - 1])
        print(('distance from surface: ' + str(depth)), file=debug_file)
    
def findInterfaceRes(pdb, interface):
    components = [[a for a in interface if a.chain in part] for part in pdb.interfaceParts]

    partResIds = [set((a.chain, a.resId) for a in part) for part in components]
    interfaceRes = [[a for a in pdb.atoms if (a.chain, a.resId) in partResIds] for part, partResIds
                    in zip(pdb.interfaceParts, partResIds)]
    return interfaceRes