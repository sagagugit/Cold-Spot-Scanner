#!/usr/bin/env python3
import os
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.ResidueDepth import get_surface, min_dist

from .. import setupPpiDb
from ..setupPpiDb import PDBS_DIR, getInterfaceAtoms
from pdbReader import PDBReader

def interaceAtomDepths(pdbName):
    filename = os.path.join(PDBS_DIR, '%s_FH.pdb' % pdbName[:4])
    pdb = PDBReader.readFile(filename, pdbsNamesToChains[pdbName[0:4]])
    interfaceAtoms = getInterfaceAtoms(cursor, pdb)

    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure(pdbName, filename)
    model = structure[0]
    chains = model.get_list()
    atoms = [chain.get_list() for chain in chains]
    interfaceAtoms = pdbReader_to_BioPyth(interfaceAtoms, Bioatoms)
    
    surface = get_surface(model)
    
    InterfaceDistances = []
    for atom in InterfaceAtoms:
        dist = min_dist(atom.get_coord(), surface)
        residue = atom.get_parent()
        chain = residue.get_parent()
        InterfaceDistance.append([chain.get_id(), residue.get_id(), atom.get_name(), dist])
    return InterfaceDistances
    
def pdbReader_to_BioPyth(interfaceAtoms, BioAtoms):
    interfaceBioAtoms = []
    for a in interfaceAtoms:
        for atom in BioAtoms:
            residue = atom.get_parent()
            if residue.get_id() == a.resId and atom.get_name() == a.symbol:
                interfaceBioAtoms.append(atom)       
    return interfaceBioAtoms