"""
Calculates salt bridges or electrostatic between two chains.

Arguments:
    saltBridges - use it to get stats for pair of 2 charges in distance of less than 4A.
    periphery - use it to look for periphery charges only

See also:
DRIEDING: A Generic Force field for molecular simulations
"""

from __future__ import print_function
import os
import numpy as np
import math
from scipy.linalg import norm
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.vectors import calc_angle

from . import DBConfig
from .kdtree import KDTree
from .ASA import ASA

VERBOSE = True

def getHbonds(pdb, pdbName):
    conn = DBConfig.get_connection()
    cursor = conn.cursor()
    cursor.execute("""select DonorChain,DonorResId,DonorICode,DonorSymbol,AccChain,AccResId,AccICode,AccSymbol from
                        Ndrieding
                        inner join donors2
                        on DonorSymbol=donors2.Symbol
                        where PDB='%s'""" % pdbName)

    hbonds = set()
    for dChain, dResId, dICode, dSymbol, aChain, aResId, aICode, aSymbol in cursor.fetchall():
        donorAtom, accAtom = None, None
        for a in pdb.atoms:
            if a.chain == dChain and a.resId == dResId and a.iCode == dICode and a.symbol == dSymbol:
                donorRes = (a.resId, a.iCode)
            elif a.chain == aChain and a.resId == aResId and a.iCode == aICode and a.symbol == aSymbol:
                accRes = (a.resId, a.iCode)
        hbonds.add((donorRes, accRes))
    return hbonds
    print(hbonds)

def fullHbonds(hbonds,charged_atom):
    charged_atomHbonds= [Hbond for Hbond in hbonds if charged_atom.resId in Hbond]
    if (charged_atom.residue in ['ASP', 'GLU']) and len(charged_atomHbonds) >= 5:
        return fullHbonds
    if charged_atom.residue == 'LYS' and len(charged_atomHbonds) >= 3:
        return fullHbonds
    if charged_atom.residue == 'ARG' and len(charged_atomHbonds) >= 6:
        return fullHbonds


def eInteraction(Qi, Qj, R):
    kcal_mol_constant = 322.0637
    return kcal_mol_constant * Qi * Qj / (R ** 2)


def assignCharge(atom, pH=7):
    # how do we assign?

    if atom.residue in ['ASP', 'GLU'] and atom.atomType == 'O' and atom.symbol != 'O':
        return -0.5  # -1
    # arg deloclalized
    # ARG - NH1 NH2 (not NE and N)
    # LYS NZ
    # HIS ND1 NE2
    if atom.residue == 'LYS' and atom.symbol == 'NZ':
        return 1.0
    posRes = ['ARG']
    if 0.1 < pH <= 6:
        posRes.append('HIS')
    if atom.residue in posRes and atom.atomType == 'N' and atom.symbol not in ['N', 'NE']:
        return 0.5  # 1

    return 0


def calcElectrostatic(pdb, interface, exclude_hbonds=False, count_cutoff_distance=4.5):
    """
    Calculated possible electro interactions, excluding 1-2 and 1-3 interactions
    (already included in angle and bond interactions
     """
    CUTOFF_DISTANCE = 4.5  # we could have 6?
    if exclude_hbonds:
        hHbonds = getHbonds(pdb, pdb.name)
    else:
        hHbonds = []

    components = []

    for part in pdb.interfaceParts:
        components.append([a for a in interface if a.chain in part and assignCharge(a) != 0])

    # fill kd tree with charged atoms from second component
    comp2 = [atom for atom in pdb.atoms if atom.chain in pdb.interfaceParts[1] and assignCharge(atom) != 0]

    ktree = KDTree.construct_from_data(comp2)
    electroStat = 0.0
    pp, pm, mm = 0, 0, 0
    for atom in components[0]:
        # interactions are not calculated between atoms bonded to each other (1,2 and 1,3 [hbonds])
        Qi = assignCharge(atom)
        if Qi == 0:
            continue

        nearAtoms = list(ktree.findByDistance(query_point=atom.coord, distance=CUTOFF_DISTANCE ** 2))
        contact = [con for con in nearAtoms if not (((con, atom) in hHbonds) or ((atom, con) in hHbonds))]
        for con in contact:
            Qj = assignCharge(con)
            if Qj == 0:
                continue
            R = math.sqrt(atom.distance(con))
            electroStat += eInteraction(Qi, Qj, R)

            if R < count_cutoff_distance:
                if Qi > 0 and Qj > 0:
                    pp += 1
                elif Qi < 0 and Qj < 0:
                    mm += 1
                else:
                    pm += 1

    return electroStat, pp, mm, pm


def is_hydrophilic(atom):
    """
    Checks whether an atom belongs to hydrophilic residue
    :param atom: atom
    :return: True if the atom belongs to hydrophobic residue, otherwise false
    """

    hydrophilic_atoms = ['H', 'N', 'S', 'O']
    hydrophilic_residues = ['GLU', 'ASP', 'ASN', 'QLN', 'HIS', 'GLN', 'SER', 'THR', 'ARG', 'LYS']
    
    if atom.atomType in hydrophilic_atoms:
        return is_hydrophilic
    if atom.residue in hydrophilic_residues:
        return is_hydrophilic
    if atom.symbol in ['C', 'CA']:
        return is_hydrophilic
    if atom.residue == 'TRP' and atom.symbol in ['CE2', 'CD1']:
        return is_hydrophilic    
    if atom.residue == 'PRO' and atom.symbol == 'CD':
        return is_hydrophilic
    if atom.residue == 'TYR' and atom.symbol == 'CZ':
        return is_hydrophilic
    if atom.residue == 'MET' and atom.symbol != 'CE':
        return is_hydrophilic


def find_charged_pi(pdb, interface, depthDistances):
    
    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure(pdb.name, pdb.file)
    model = structure[0]
    
    HYDROPHOBIC_CHARGED_CUTOFF_DISTANCE = 4.5  # we could have 6?

    Hbonds = getHbonds(pdb, pdb.name)


    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure(pdb.name, pdb.file)
    model = structure[0]

    charged_pi_interactions = []
    non_charged_pi_interactions = []

    for part in pdb.interfaceParts:
        charged_atoms = [a for a in interface if a.residue in ["LYS", "ARG"] and assignCharge(a) != 0 and get_residue_depth(a, depthDistances) >= 4.0]
        if not any(charged_atoms):
            continue

            electro_kdtree = KDTree.construct_from_data(charged_atoms)
            hydrophobic_partners = [a for a in interface if a.residue in ["TYR", "TRP", "PHE"] and not is_hydrophilic(a)]
            for atom in hydrophobic_partners:
                nearAtoms = list(electro_kdtree.findByDistance(query_point=atom.coord, distance=HYDROPHOBIC_CHARGED_CUTOFF_DISTANCE ** 2))
                for con in nearAtoms:
                    R = math.sqrt(atom.distance(con))
                    depth = get_residue_depth(con, depthDistances)   
                    if R < HYDROPHOBIC_CHARGED_CUTOFF_DISTANCE:
                        if ((con.resId, con.iCode), (atom.resId, atom.iCode)) in Hbonds or ((atom.resId, atom.iCode), (con.resId, con.iCode)) in Hbonds:
                            continue
                
                    if is_charged_pi(con, atom, model, pdb):
                        charged_pi_interactions.append((con, atom))
                    else:
                        non_charged_pi_interactions.append((con, atom))

    return charged_pi_interactions

def is_charged_pi(con, atom, model, pdb):

    atom_res = model[atom.chain][atom.resId]
    if con.residue == "LYS":
        con_atom = model[con.chain][con.resId][con.symbol].get_coord()
    if con.residue == "ARG":
        con_atom = model[con.chain][con.resId]["CZ"].get_coord()
    
    if atom_res.get_resname() in ["PHE", "TYR"]:
        
        CG_atom = atom_res["CG"].get_coord()
        CZ_atom = atom_res["CZ"].get_coord()
        CD1_atom = atom_res["CD1"].get_coord()
        
        midpoint = ( CG_atom + CZ_atom ) / 2
        
        vect1 =  CZ_atom - midpoint
        vect2 =  CD1_atom - midpoint
        
    if atom_res.get_resname() == "TRP":
        
        CE2_atom = atom_res["CE2"].get_coord()
        CE3_atom = atom_res["CE3"].get_coord()
        CH2_atom = atom_res["CH2"].get_coord()
            
        midpoint = ( CE3_atom + CH2_atom ) / 2
        
        vect1 = CH2_atom - midpoint
        vect2 = CE3_atom - midpoint
        
    norm_vect = np.cross(vect1, vect2)
    norm_vect_abs = norm(norm_vect)
    #the cross product of two vectors is a vector that is 90 degrees to the plane formed by the two vectors
    
    vect3 = con_atom - midpoint
    vect3_abs = norm(vect3)

    angle = np.rad2deg(np.arccos(np.dot(norm_vect, vect3)/(norm_vect_abs * vect3_abs)))
    #uses definition of dot product ( a dot b = a * b * cos(theta) ) where theta is the angles between vectors a and b
    
    if vect3_abs > 6:
        return False

    if (abs(angle) < 45) or (abs(angle) > (45 + 90)):
        return True
    else:
        return False



def anion_aromatic(pdb, interface, depthDistances):
    
    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure(pdb.name, pdb.file)
    model = structure[0]
    
    anion_aromatic_distance = 6  # we could have 6?

    Hbonds = getHbonds(pdb, pdb.name)


    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure(pdb.name, pdb.file)
    model = structure[0]

    anion_aromatic_interactions = []
    

    anion_aromatic_atoms = [a for a in interface if a.residue in ["ASP", "GLU"] and assignCharge(a) != 0]
    
    if not any(anion_aromatic_atoms):
        return anion_aromatic_interactions

    electro_kdtree = KDTree.construct_from_data(anion_aromatic_atoms)
    aromatic_partners = [a for a in interface if a.residue in ["TYR", "TRP", "PHE"] and not is_hydrophilic(a)]

    for atom in aromatic_partners:
        nearAtoms = list(electro_kdtree.findByDistance(query_point=atom.coord, distance=anion_aromatic_distance ** 2))
        for con in nearAtoms:
            R = math.sqrt(atom.distance(con))
            depth = get_residue_depth(con, depthDistances)
            if R < anion_aromatic_distance and not (atom.chain == con.chain and atom.resId == con.resId and atom.iCode == con.iCode):
                if ((con.resId, con.iCode), (atom.resId, atom.iCode)) in Hbonds or ((atom.resId, atom.iCode), (con.resId, con.iCode)) in Hbonds:
                    continue
                if is_anion_aromatic_possible(con, atom, model, pdb):
                    anion_aromatic_interactions.append((con, atom))

    return anion_aromatic_interactions



def is_anion_aromatic_possible(neg, atom, model, pdb):
    global vect11, vect22, midpoint, neg_atom

    atom_res = model[atom.chain][atom.resId]
    try:
        if neg.residue == "ASP":
            neg_atom = model[neg.chain][neg.resId]["OD1" or "OD2"].get_coord()
        if neg.residue == "GLU":
            neg_atom = model[neg.chain][neg.resId]["OE1" or "OE2"].get_coord()
    except:
        pass
    
    if atom_res.get_resname() in ["PHE", "TYR"]:
        
        CG_atom = atom_res["CG"].get_coord()
        CZ_atom = atom_res["CZ"].get_coord()
        CD1_atom = atom_res["CD1"].get_coord()
        
        midpoint = ( CG_atom + CZ_atom ) / 2
        
        vect11 =  CZ_atom - midpoint
        vect22 =  CD1_atom - midpoint
        
    if atom_res.get_resname() == "TRP":
        
        CE2_atom = atom_res["CE2"].get_coord()
        CE3_atom = atom_res["CE3"].get_coord()
        CH2_atom = atom_res["CH2"].get_coord()
            
        midpoint = ( CE3_atom + CH2_atom ) / 2
        
        vect11 = CH2_atom - midpoint
        vect22 = CE3_atom - midpoint
        
    norm_vect = np.cross(vect11, vect22)
    norm_vect_abs = norm(norm_vect)
    #the cross product of two vectors is a vector that is 90 degrees to the plane formed by the two vectors
    
    vect3 = neg_atom - midpoint
    vect3_abs = norm(vect3)

    angle = np.rad2deg(np.arccos(np.dot(norm_vect, vect3)/(norm_vect_abs * vect3_abs)))
    #uses definition of dot product ( a dot b = a * b * cos(theta) ) where theta is the angles between vectors a and b
    
    if vect3_abs > 6:
        return False

    if 80 <= (abs(angle)) <= 100:
        return True
    else:
        return False


def charged_atom_areas(atom, pdb):
    
    if atom.residue == 'ASP':
        charged = [a for a in pdb.atoms if a.symbol in ['OD1', 'OD2'] and a.resId == atom.resId]
    if atom.residue == 'GLU':
        charged = [a for a in pdb.atoms if a.symbol in ['OE1', 'OE2'] and a.resId == atom.resId]
    if atom.residue == 'LYS':
        charged = [a for a in pdb.atoms if a.symbol in ['1HZ', '2HZ', '3HZ', 'NZ'] and a.resId == atom.resId]
    if atom.residue == 'ARG':
        charged = [a for a in pdb.atoms if a.symbol in ['NH1', 'NH2', 'NE'] and a.resId == atom.resId]
    return charged

def get_residue_depth(atom, depthDistances):
    for a, dist in depthDistances:
        if a.resId == atom.resId and a.iCode == atom.iCode and a.symbol == atom.symbol:
            return dist

def calcInterElectroHydrophobic(pdb, interface, depthDistances):
    """
    Calculated possible electro interactions, excluding 1-2 and 1-3 interactions
    (already included in angle and bond interactions
     """
    HYDROPHOBIC_CHARGED_CUTOFF_DISTANCE = 4.5  # we could have 6?
    
    Hbonds = getHbonds(pdb, pdb.name)
    

    hydrophobic_charged_interactions = []
    hydroElectroOutput = pdb.getFile('.hydrophobicElectro.txt') if VERBOSE else None

    for part in pdb.interfaceParts:
        charged_atoms = [a for a in interface if assignCharge(a) != 0 and a.chain in part and get_residue_depth(a, depthDistances) >= 4.0]
        if not any(charged_atoms):
            continue
        electro_kdtree = KDTree.construct_from_data(charged_atoms)
        other_parts = ''.join([partb for partb in pdb.interfaceParts if partb != part])
        hydrophobic_partners = [a for a in interface if a.chain in other_parts and not is_hydrophilic(a)]
        for atom in hydrophobic_partners:
            nearAtoms = list(electro_kdtree.findByDistance(query_point=atom.coord, distance=HYDROPHOBIC_CHARGED_CUTOFF_DISTANCE ** 2))
            for con in nearAtoms:
                Qi = assignCharge(con)
                R = math.sqrt(atom.distance(con))
                depth = get_residue_depth(con, depthDistances)
                if R < HYDROPHOBIC_CHARGED_CUTOFF_DISTANCE:
                    if ((con.resId, con.iCode), (atom.resId, atom.iCode)) in Hbonds or ((atom.resId, atom.iCode), (con.resId, con.iCode)) in Hbonds:
                        continue
                if Qi > 0:
                    hydrophobic_charged_interactions.append((con, atom))
                elif Qi < 0:
                    hydrophobic_charged_interactions.append((con, atom))
                if hydroElectroOutput:
                    print(','.join((con.chain, str(con.resId), con.residue, atom.chain, str(atom.resId),
                                    atom.residue, '%.3f' % R)), file=hydroElectroOutput)
    if hydroElectroOutput:
        hydroElectroOutput.close()
    return hydrophobic_charged_interactions    

def calcIntraElectroHydrophobic(pdb, interface, depthDistances):
    """
    Calculated possible electro interactions, excluding 1-2 and 1-3 interactions
    (already included in angle and bond interactions
     """
    HYDROPHOBIC_CHARGED_CUTOFF_DISTANCE = 4.5  # we could have 6?
    
    Hbonds = getHbonds(pdb, pdb.name)
    

    hydrophobic_charged_interactions2 = []
    hydroElectroOutput = pdb.getFile('.hydrophobicElectro.txt') if VERBOSE else None

    for part in pdb.interfaceParts:
        charged_atoms = [a for a in interface if assignCharge(a) != 0 and a.chain in part and
                         get_residue_depth(a, depthDistances) >= 4.0]
        if not any(charged_atoms):
            continue
        electro_kdtree = KDTree.construct_from_data(charged_atoms)
        hydrophobic_partners = [a for a in interface if a.chain in part and not is_hydrophilic(a)]
        for atom in hydrophobic_partners:
            nearAtoms = list(electro_kdtree.findByDistance(query_point=atom.coord, distance=HYDROPHOBIC_CHARGED_CUTOFF_DISTANCE ** 2))
            for con in nearAtoms:
                depth = get_residue_depth(con, depthDistances)
                Qi = assignCharge(con)
                R = math.sqrt(atom.distance(con))
                if R < HYDROPHOBIC_CHARGED_CUTOFF_DISTANCE and not (atom.resId == con.resId and atom.iCode == con.iCode):
                    if (((con.resId, con.iCode), (atom.resId, atom.iCode)) in Hbonds or ((atom.resId, atom.iCode), (con.resId, con.iCode)) in Hbonds):
                        continue
                    if Qi > 0:
                        hydrophobic_charged_interactions2.append((con, atom))
                    elif Qi < 0:
                        hydrophobic_charged_interactions2.append((con, atom))
                    if hydroElectroOutput:
                        print(','.join((con.chain, str(con.resId), con.residue, atom.chain, str(atom.resId),
                                       atom.residue, '%.3f' % R)), file=hydroElectroOutput)
    if hydroElectroOutput:
        hydroElectroOutput.close()
    return hydrophobic_charged_interactions2    

def find_charged_Negative(pdb, interface):
    
    Charge_CUTOFF_DISTANCE = 3.5 # we could have 6?


    charge_interactions = []
    positive_atoms = [a for a in interface if a.residue in ['LYS', 'ARG'] and assignCharge(a) != 0  and is_positive(a)]
    if not any(positive_atoms):
        return charge_interactions

    electro_kdtree = KDTree.construct_from_data(positive_atoms)
    negative_partners = [a for a in interface if a.residue in ['ASP', 'GLU'] and assignCharge(a) != 0  and is_negative(a)]

    for atom in negative_partners:
        nearAtoms = list(electro_kdtree.findByDistance(query_point=atom.coord, distance=Charge_CUTOFF_DISTANCE ** 2))
        for con in nearAtoms:
            Qi = assignCharge(con)
            R = math.sqrt(atom.distance(con))
            if R < Charge_CUTOFF_DISTANCE:
                if Qi > 0:
                    charge_interactions.append((con, atom))
                elif Qi < 0:
                    charge_interactions.append((con, atom))
                    
                
    return charge_interactions



def is_positive(atom):
    if atom.residue == 'ARG' and atom.symbol in ['NH1','NH2']:
        return is_positive    
    if atom.residue == 'LYS' and atom.symbol == 'NZ':
        return is_positive

    
def is_negative(atom):
    if atom.residue == 'ASP' and atom.symbol in ['OD1', 'OD2']:
        return is_negative
    if atom.residue == 'GLU' and atom.symbol in ['OE1', 'OE2']:
        return is_negative

def identifying_interface_ligands(pdb, interface):
    """Check whether the ligand are in a distance of 3A from the interface atoms"""

    ligans = []

    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure(pdb.name, pdb.file)
    model = structure[0]
    

    all_interface_atoms = [a for a in interface]
    
    for con in all_interface_atoms:
        if is_ligand_interface(con, pdb):
            ligans.append(con)

    return ligans
                
    

def is_ligand_interface(con, pdb):

    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure(pdb.name, pdb.file)
    model = structure[0]

    list_c = []

    for con in structure.get_atoms():
        list_c.append(con.get_coord())

    ligand_distances = []


    for atoms in structure.get_atoms():
        tags = atoms.get_full_id()
        if tags[3][0] != " ":
            if tags[3][0] != 'W':
                for x in list_c:
                    if (atoms.get_coord() - x).any() < 4:
                        return True
                    else:
                        False

def identifying_same_charge(pdb, interface):
    """identifying whether there are same charge atoms 
    near by to each other in the interface"""

    same_charge_cut_off_distance = 4.5

    Hbonds = getHbonds(pdb, pdb.name)


    same_charge_interactions = []
    hydroElectroOutput = pdb.getFile('.hydrophobicElectro.txt') if VERBOSE else None

    for part in pdb.interfaceParts:
        first_positive_atom = [a for a in interface if a.chain in part and is_positive(a)]
        if not any(first_positive_atom):
            continue
        electro_kdtree = KDTree.construct_from_data(first_positive_atom)
        other_parts = ''.join([partb for partb in pdb.interfaceParts if partb != part])
        second_positive_atom =[a for a in interface if a.chain in other_parts and is_positive(a)]
        for atom in second_positive_atom:
            nearAtoms = list(electro_kdtree.findByDistance(query_point=atom.coord, distance=same_charge_cut_off_distance ** 2))
            for con in nearAtoms:
                R = math.sqrt(atom.distance(con))
                if R < same_charge_cut_off_distance:
                    if ((con.resId, con.iCode), (atom.resId, atom.iCode)) in Hbonds or ((atom.resId, atom.iCode), (con.resId, con.iCode)) in Hbonds:
                        continue
                    same_charge_interactions.append((con, atom))


    for part in pdb.interfaceParts:
        first_negative_atom = [a for a in interface if a.chain in part and is_negative(a)]
        if not any(first_negative_atom):
            continue
        electro_kdtree = KDTree.construct_from_data(first_negative_atom)
        other_parts = ''.join([partb for partb in pdb.interfaceParts if partb != part])
        second_negative_atom =[a for a in interface if a.chain in other_parts and is_negative(a)]
        for atom in second_negative_atom:
            nearAtoms = list(electro_kdtree.findByDistance(query_point=atom.coord, distance=same_charge_cut_off_distance ** 2))
            for con in nearAtoms:
                R = math.sqrt(atom.distance(con))
                if R < same_charge_cut_off_distance:
                    if ((con.resId, con.iCode), (atom.resId, atom.iCode)) in Hbonds or ((atom.resId, atom.iCode), (con.resId, con.iCode)) in Hbonds:
                        continue
                    same_charge_interactions.append((con, atom))

    final_same_charge_interactions = list(set(same_charge_interactions))

    return final_same_charge_interactions










