#!/usr/bin/env python
from __future__ import print_function
import argparse
import os
import subprocess
import sys
import pkg_resources
import requests
import numpy as np
from Bio.PDB.PDBParser import PDBParser
from pyPPI.ResidueDepthcopy import get_surface, min_dist
from pyPPI import DBConfig
import pyPPI.surfaceComplementarity.VDW as VDW
import pyPPI.surfaceComplementarity.interfaceDepth as Periphery
from pyPPI.ASA import ASA
from pyPPI.hbonds import hbonds
from pyPPI.kdtree import KDTree
import pyPPI.pdbReader as pdbReader
from pyPPI.pdbReader import PDBReader
import pyPPI.electrostat as electrostat
import pyPPI.cavities as cavities
import pyPPI.Solvent_inaccessible_cavities as si_cavities
import pyPPI.heteroatom as heteroatom
from itertools import groupby
from pprint import pprint
from collections import defaultdict
from Bio import PDB
from Bio import PDB
from sklearn.cluster import DBSCAN
import numpy as np
from collections import Counter
import os
from Bio.PDB.PDBParser import PDBParser
from pyPPI.PDBIO import Select, PDBIO
import csv
import math
import logging
import stat
import pandas as pd
"""
Distance in angtroms between the chains that is relevant for defining the interface
"""
INTERFACE_DISTANCE = 4
WORKING_DIRECTORY = './'
PDBS_DIR = "./pdbs/"
RESULTS_DIR = "./results/"
PARAMETERS_DIR = "./parameters/"
COLD_SPOTS_DIR = "./cold_spots/"

_remediator = pkg_resources.resource_filename('pyPPI', '/'.join(['molprobity', 'remediator.pl']))
_reduce_path = pkg_resources.resource_filename('pyPPI', '/'.join(['molprobity', 'reduce']))


def download_PDB(pdb):
    """
    Downloads a PDB from protein data base
    :param pdb: pdb identifier
    """
    url = 'http://www.rcsb.org/pdb/files/{0}.pdb'.format(pdb)
    print('downloading %s (%s)' % (pdb, url))

    req = requests.get(url)
    with get_file(pdb) as newPDB:
        print(req.text, file=newPDB)


def get_file(name):
    """
    Get file for write in the PDBS_DIR
    :param name:
    :return:
    """
    global PDBS_DIR
    return open(os.path.join(name + ".pdb"), "w")


def download_DB(pdbList):
    """
    Downloads PDB and add hydrogens using molprobity
    :param pdbList: list of pdbs to download
    """
    print("Downloading pdbs according to list")
    for pdb in pdbList:
        # don't download twice the same PDB
        if os.path.exists(os.path.join(pdb + "_FH.pdb")): continue

        # in case the PDB is already in the directory
        if not os.path.exists(os.path.join(pdb + ".pdb")):
            download_PDB(pdb)

        molprobity(pdb)
    print("Finished downloading pdbs")


def molprobity(pdb_name):
    """
    runs molprobility on a input protein
    :param pdb_name: name of the PDB file
    :return:
    """
    global MOLPROBITY_DIR, PDBS_DIR
    if os.path.exists(os.path.join(pdb_name + "_FH.pdb")):
        return True  # already exist
    print('Starting molprobity %s' % pdb_name)
    subprocess.check_output('perl ' + _remediator + ' ' + os.path.join(
                                                                                                         pdb_name + ".pdb") + ' > a%s' % pdb_name,
                            shell=True)
    try:
        subprocess.check_output('./pyPPI/molprobity/reduce' + ' a%s> b%s' % (pdb_name, pdb_name) , shell=True)
    except:
        print('error prasing PDB %s' % pdb_name)
        pass  # yakky kaky, but reduce returns 1 exit
    subprocess.check_output(
        'perl ' + _remediator +' b%s -oldout> '% pdb_name + os.path.join(PDBS_DIR, pdb_name + "_FH.pdb"),
        shell=True)
    # delete the PDB file - we will work with a file with hydrogens added (_FH create above)
    #os.remove(os.path.join(pdb_name + ".pdb"))


def buildASAperAtomForComplex(pdb, result):
    asaCalc = ASA(pdb)
    asaCalc.execute()
    for atom, asa in asaCalc.interPerAtom.items():
        # complex inter
        res = [pdb.name, atom.chain, atom.residue, atom.resId, atom.iCode, atom.symbol, atom.atomType, asa, atom.tempFactor, 0]
        print(','.join([str(a) for a in res]), file=result)
        # complex intra (separated)
        asa = asaCalc.diffASAperAtom[atom] + asa
        res = [pdb.name, atom.chain, atom.residue, atom.resId, atom.iCode, atom.symbol, atom.atomType, asa, atom.tempFactor, 1]
        print(','.join([str(a) for a in res]), file=result)

def calcInterfaceDist(pdb, result):
    global INTERFACE_DISTANCE
    partA = [a for a in pdb.atoms if a.chain in pdb.interfaceParts[0]]
    partB = [a for a in pdb.atoms if a.chain in pdb.interfaceParts[1]]
    if len(partA) == 0 or len(partB) == 0:
        print('WARNING: %s doesnt have atoms in one its chains' % pdb.name)
        return
    aTree = KDTree.construct_from_data(partA[:])
    bTree = KDTree.construct_from_data(partB[:])
    complexChains = ':'.join(pdb.interfaceParts)
    for part, tree in [(partA, bTree), (partB, aTree)]:
        for atom in part:
            near, dist = tree.findNearest(query_point=atom.coord, num=1)
            if dist < INTERFACE_DISTANCE:
                print(','.join([pdb.name, complexChains, atom.chain, str(atom.resId), atom.iCode, atom.symbol, atom.atomType, str(dist)]), file=result)

def createInterfaceCSV(pdbsToAnalyze):
    global PDBS_DIR, RESULTS_DIR, PARAMETERS_DIR
    if all(os.path.exists(os.path.join(RESULTS_DIR, PARAMETERS_DIR, resFile)) for resFile in ['PerAtomASA.csv', 'PerAtomASA.csv']):
        print('Data already exist in result directory.')
        return
    
    with open(os.path.join(RESULTS_DIR, PARAMETERS_DIR, 'PerAtomASA.csv'), 'w') as asaPerAtom:
        with open(os.path.join(RESULTS_DIR, PARAMETERS_DIR, 'PerAtomDistance.csv'), 'w') as distancePerAtom:
            pdbs = os.listdir(PDBS_DIR)
            failedPDBs = []
            pdbsNamesToChains = dict((p[0], p[1].split(':') if len(p) > 1 else None) for p in pdbsToAnalyze)
            for pdbName in pdbs:
                if pdbName[0:4] not in pdbsNamesToChains: continue
                pdb = PDBReader.readFile(os.path.join(PDBS_DIR, pdbName), pdbsNamesToChains[pdbName[0:4]])
                try:
                    print('Writing ASA for %s' % pdb.name)
                    buildASAperAtomForComplex(pdb, asaPerAtom)
                    print('Writing distance for %s' % pdb.name)
                    calcInterfaceDist(pdb, distancePerAtom)
                except IndexError:
                    failedPDBs.append(pdb.name)

    print('Finished')
    if len(failedPDBs) > 0:
        print('Failed to process:', ','.join(failedPDBs))

class StdevFunc:
    def __init__(self):
        self.M = 0.0
        self.S = 0.0
        self.k = 1

    def step(self, value):
        if value is None:
            return
        tM = self.M
        self.M += (value - tM) / self.k
        self.S += (value - tM) * (value - self.M)
        self.k += 1

    def finalize(self):
        if self.k < 3:
            return None
        return math.sqrt(self.S / (self.k-2))



def getInterfaceAtoms(cur, pdb):
    cur.execute('''
    select Chain,ResId,iCode,Symbol from NinterfaceAtoms
    where PDB='%s'
    ''' % pdb.name)
    interfaceAtoms = []
    for chain, resid, icode, symbol in cur.fetchall():
        interfaceAtoms.append(next(a for a in pdb.atoms if a.chain == chain and a.resId == resid and a.iCode == icode and a.symbol == symbol))
    return interfaceAtoms

class ChainSelect(Select):
    def __init__(self, chain):
        self.chain = chain

    def accept_chain(self, chain):
        if chain.get_id() in self.chain:
            return 1
        else:          
            return 0
class NonHetSelect(Select):
    def accept_residue(self, residue):
        return 1 if residue.id[0] == " " else 0


def getInterfaceHeteroatoms(cur, pdb):
    cur.execute('''
    select Chain,ResId,iCode,Symbol from NinterfaceAtoms
    where PDB='%s'
    ''' % pdb.name)
    interfaceheteroAtoms = []
    for c in cur.fetchall():
        interfaceheteroAtoms.append(next(a for a in pdb.hetAtms))
    return interfaceheteroAtoms

def give_permission_msms():
    st = os.stat('./pyPPI/msms/msms')
    return os.chmod('./pyPPI/msms/msms', st.st_mode | stat.S_IEXEC)

def give_permission_reduce():
    st = os.stat('./pyPPI/molprobity/reduce')
    return os.chmod('./pyPPI/molprobity/reduce', st.st_mode | stat.S_IEXEC)


def fillInterfacePeriphrial(pdbsToAnalyze):
    global PDBS_DIR, RESULTS_DIR, PARAMETERS_DIR

    if os.path.exists(os.path.join(RESULTS_DIR, PARAMETERS_DIR, 'interfacePeriphrial.csv')):
        print('Data already exist in result directory for interface periphery.')
        return

    pdbsNamesToChains = dict((p[0], p[1].split(':') if len(p) > 1 else None) for p in pdbsToAnalyze)

    with open(os.path.join(RESULTS_DIR, PARAMETERS_DIR, 'interfacePeriphrial.csv'), 'w') as interfacePeriphrial:
    #    print('PDB,Chain,ResId,Symbol,Peripherial,PropPeri', file=interfacePeriphrial)
        for pdbName, chains in pdbsNamesToChains.items():
            print('Calculating peripheral table for %s ' % pdbName)
            pdb_path = os.path.join(PDBS_DIR, '%s_FH.pdb' % pdbName)
            depthL, peripherialL = Periphery.calc_peripheral_PDB(pdb_path, chains)
            for atom, peri, propPeri in peripherialL:
                print(','.join([pdbName, atom.chain, str(atom.resId), atom.symbol, str(peri), str(propPeri)]),
                      file=interfacePeriphrial)

    conn = DBConfig.get_connection()
    cursor = conn.cursor()
    with open(os.path.join(RESULTS_DIR,PARAMETERS_DIR, 'interfacePeriphrial.csv'), 'r') as csv_data2:
        csvfile2 =csv.reader(csv_data2, delimiter=',')
        for row2 in csvfile2:
            value2 = (row2[0], row2[1], row2[2], row2[3], row2[4], row2[5])            
            cursor.execute("INSERT INTO `interfacePeriphrial` (`PDB`,`Chain`,`ResId`,`Symbol`,`Peri`,`PropPeri`) values (?,?,?,?,?,?);", value2)
    conn.commit()
    conn.close()


def residue_depth(pdb):    
    pdbName = pdb.name
    filename = pdb.file
    ReaderAtomsInput = pdb.atoms
    
    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure(pdbName, filename)
    model = structure[0]
    
    BioAtoms = model.get_atoms()
    
    surface = get_surface(model, probe_radius=1.5)
    BioDepthDistances = []
    for atom in BioAtoms:
        dist = min_dist(atom.get_coord(), surface)
        BioDepthDistances.append([atom, dist])
    pdbReaderDistances = BioPyth_to_pdbReader(BioDepthDistances, ReaderAtomsInput, pdb)
    return pdbReaderDistances

def BioPyth_to_pdbReader(BioAtomsInput, ReaderAtomsInput, pdb):
    ReaderAtomsOutput = []
    for atom, dist in BioAtomsInput:
        for a in ReaderAtomsInput:
            res = atom.get_parent()
            chain = res.get_parent()
            if chain.get_id() == a.chain and res.get_id() == (' ', a.resId, a.iCode) and atom.get_name() == a.symbol:
                ReaderAtomsOutput.append([a, dist])
                break
    return ReaderAtomsOutput


def get_electrohydro_dicts(pdbsNamesToChains, cursor):
    pdb_to_inter_interactions = {}
    pdb_to_intra_interactions = {}
    for pdbName, chains in pdbsNamesToChains.items():
        pdb_path = os.path.join(PDBS_DIR, '%s_FH.pdb' % pdbName)
        pdb = PDBReader.readFile(pdb_path, chains)
        interfaceAtoms = getInterfaceAtoms(cursor, pdb)
        depthDistances = residue_depth(pdb)
        pdb_to_inter_interactions[pdbName[:4]] = electrostat.calcInterElectroHydrophobic(pdb, interfaceAtoms, depthDistances)
        pdb_to_intra_interactions[pdbName[:4]] = electrostat.calcIntraElectroHydrophobic(pdb, interfaceAtoms, depthDistances)
    return pdb_to_inter_interactions, pdb_to_intra_interactions

def get_pi_dicts(pdbsNamesToChains, cursor):
    pdb_to_pi_interactions = {}
    for pdbName, chains in pdbsNamesToChains.items():
        pdb_path = os.path.join(PDBS_DIR, '%s_FH.pdb' % pdbName)
        pdb = PDBReader.readFile(pdb_path, chains)
        interfaceAtoms = getInterfaceAtoms(cursor, pdb)
        depthDistances = residue_depth(pdb)
        pdb_to_pi_interactions[pdbName[:4]] = electrostat.find_charged_pi(pdb, interfaceAtoms, depthDistances)
    return pdb_to_pi_interactions

def get_pi_dicts_subtract(pdbsNamesToChains, cursor):
    pdb_to_pi_interactions = {}
    positive_charge_aminos = [] 
    for pdbName, chains in pdbsNamesToChains.items():
        pdb_path = os.path.join(PDBS_DIR, '%s_FH.pdb' % pdbName)
        pdb = PDBReader.readFile(pdb_path, chains)
        interfaceAtoms = getInterfaceAtoms(cursor, pdb)
        depthDistances = residue_depth(pdb)
        pdb_to_pi_interactions[pdbName[:4]] = electrostat.find_charged_pi(pdb, interfaceAtoms, depthDistances)
    for pdbName, chains in pdbsNamesToChains.items():
        if pdbName[0:4] not in pdbsNamesToChains: continue
        inter_pi_interactions = pdb_to_pi_interactions[pdbName[:4]]
        for charged_atom,  hydrophobic_atom in inter_pi_interactions:
            positive_charge_aminos.append([pdbName, charged_atom.chain, str(charged_atom.resId), charged_atom.residue, charged_atom.symbol])

    return positive_charge_aminos

def get_anionic_dicts(pdbsNamesToChains, cursor):
    pdb_to_anionic_interactions = {}
    for pdbName, chains in pdbsNamesToChains.items():
        pdb_path = os.path.join(PDBS_DIR, '%s_FH.pdb' % pdbName)
        pdb = PDBReader.readFile(pdb_path, chains)
        interfaceAtoms = getInterfaceAtoms(cursor, pdb)
        depthDistances = residue_depth(pdb)
        pdb_to_anionic_interactions[pdbName[:4]] = electrostat.anion_aromatic(pdb, interfaceAtoms, depthDistances)
    return pdb_to_anionic_interactions



def get_anion_aromatic_subtract(pdbsNamesToChains, cursor):
    pdb_to_anion_aromatic_interactions = {}
    negatively_charge_aminos = [] 
    for pdbName, chains in pdbsNamesToChains.items():
        pdb_path = os.path.join(PDBS_DIR, '%s_FH.pdb' % pdbName)
        pdb = PDBReader.readFile(pdb_path, chains)
        interfaceAtoms = getInterfaceAtoms(cursor, pdb)
        depthDistances = residue_depth(pdb)
        pdb_to_anion_aromatic_interactions[pdbName[:4]] = electrostat.find_charged_pi(pdb, interfaceAtoms, depthDistances)
    for pdbName, chains in pdbsNamesToChains.items():
        if pdbName[0:4] not in pdbsNamesToChains: continue
        inter_anion_interactions = pdb_to_anion_aromatic_interactions[pdbName[:4]]
        for charged_atom,  hydrophobic_atom in inter_anion_interactions:
            negatively_charge_aminos.append([pdbName, charged_atom.chain, str(charged_atom.resId), charged_atom.residue, charged_atom.symbol])

    return negatively_charge_aminos


def get_unfavorable_interactions(pdbsNamesToChains, cursor):
    pdb_to_unfavorable_interactions = {}
    for pdbName, chains in pdbsNamesToChains.items():
        pdb_path = os.path.join(PDBS_DIR, '%s_FH.pdb' % pdbName)
        pdb = PDBReader.readFile(pdb_path, chains)
        interfaceAtoms = getInterfaceAtoms(cursor, pdb)
        depthDistances = residue_depth(pdb)
        inter_hydrophobic = electrostat.calcInterElectroHydrophobic(pdb, interfaceAtoms, depthDistances)
        intra_hydrophobic = electrostat.calcIntraElectroHydrophobic(pdb, interfaceAtoms, depthDistances)
        cationpiss= electrostat.find_charged_pi(pdb, interfaceAtoms, depthDistances)
        ani_aro = electrostat.anion_aromatic(pdb, interfaceAtoms, depthDistances)
        total_hydropho = inter_hydrophobic + intra_hydrophobic
        unfavorable_catt_2 = set(total_hydropho) - set(cationpiss)
        unfavorable_catt = set(unfavorable_catt_2) -set(ani_aro)
        pdb_to_unfavorable_interactions[pdbName[:4]] = unfavorable_catt
    return pdb_to_unfavorable_interactions



def charge_positive_negative_subbstract(pdbsNamesToChains, cursor):
    pdb_to_positive_negative = {}
    positive_residue_and_chain = []
    for pdbName, chains in pdbsNamesToChains.items():
        pdb_path = os.path.join(PDBS_DIR, '%s_FH.pdb' % pdbName)
        pdb = PDBReader.readFile(pdb_path, chains)
        interfaceAtoms = getInterfaceAtoms(cursor, pdb)
        pdb_to_positive_negative[pdbName[:4]] = electrostat.find_charged_Negative(pdb, interfaceAtoms)
    for pdbName, chains in pdbsNamesToChains.items():
        if pdbName[0:4] not in pdbsNamesToChains: continue
        positive_negative = pdb_to_positive_negative[pdbName[:4]]
        for charged_atom,  negatived_atom in positive_negative:
            pos=(pdbName, charged_atom.chain, str(charged_atom.resId), charged_atom.residue, charged_atom.symbol)
            positive_residue_and_chain.append(pos)
        for charged_atom,  negative_atom in positive_negative:
            neg=(pdbName, negative_atom.chain, str(negative_atom.resId), negative_atom.residue, negative_atom.symbol)
            positive_residue_and_chain.append(neg)

    return positive_residue_and_chain


def same_charge_interface(pdbsNamesToChains, cursor):
    pdb_to_same_charge = {}
    for pdbName, chains in pdbsNamesToChains.items():
        pdb_path = os.path.join(PDBS_DIR, '%s_FH.pdb' % pdbName)
        pdb = PDBReader.readFile(pdb_path, chains)
        interfaceAtoms = getInterfaceAtoms(cursor, pdb)
        pdb_to_same_charge[pdbName[:4]]= electrostat.identifying_same_charge(pdb, interfaceAtoms)

    return pdb_to_same_charge

def inserting_the_data(pdbsToAnalyzeWithChains):
    conn = DBConfig.get_connection()
    conn.create_aggregate("stdev", 1, StdevFunc)
    cursor = conn.cursor()

    with open(os.path.join(RESULTS_DIR,PARAMETERS_DIR, 'PerAtomASA.csv'), 'r') as csv_data1:
        csvfile1 =csv.reader(csv_data1, delimiter=',')
        for row1 in csvfile1:
            value1 = (row1[0], row1[1], row1[2], row1[3], row1[4], row1[5], row1[6], row1[7], row1[8], row1[9])        
            cursor.execute("INSERT INTO `perAtomASA` (`PDB`,`Chain`,`Residue`,`ResId`,`iCode`,`Symbol`,`Atom`,`ASA`,`Bfactor`,`Seperated`) values (?,?,?,?,?,?,?,?,?,?);", value1)
    with open(os.path.join(RESULTS_DIR,PARAMETERS_DIR, 'PerAtomDistance.csv'), 'r') as csv_data:
        csvfile =csv.reader(csv_data, delimiter=',')
        for row in csvfile:
            value = (row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7])
            cursor.execute("INSERT INTO `interfaceDist` (`PDB`,`Chains`,`Chain`,`ResId`,`iCode`,`Symbol`,`Atom`,`MinDist`) values (?,?,?,?,?,?,?,?);", value)
    cursor.execute("INSERT INTO `NinterfaceAtoms` (`PDB`, `Chain`, `Residue`, `ResId`, `iCode`, `Symbol`, `atom`, `diffASA`) SELECT `PDB`,`Chain`,`Residue`,`ResId`,`iCode`,`Symbol`,`Atom`, max(ASA)-min(ASA) FROM `perAtomASA` GROUP BY PDB,Chain,Residue,ResId,iCode,Symbol,Atom  HAVING stdev(ASA)>0;")
    cursor.execute("INSERT OR IGNORE INTO `NinterfaceAtoms` (`PDB`,`Chain`,`Residue`,`ResId`,`iCode`,`Symbol`,`atom`) SELECT asa.PDB, asa.Chain, asa.Residue, asa.ResId, asa.iCode, asa.Symbol, dist.Atom from `interfaceDist` `dist` inner join perAtomASA asa on dist.PDB=asa.PDB and dist.Chain=asa.Chain and dist.ResId=asa.ResId and dist.iCode=asa.iCode and dist.Symbol=asa.Symbol and Seperated=0;")



    conn.commit()

    pdbsNamesToChains = dict((p[0], p[1].split(':') if len(p) > 1 else None) for p in pdbsToAnalyzeWithChains)
    dataToInsert = []
    for pdbName, chains in pdbsNamesToChains.items():
        pdb = PDBReader.readFile(os.path.join(PDBS_DIR, '%s_FH.pdb' % pdbName), pdbsNamesToChains[pdbName[0:4]])
        if chains is None:
            compunds = pdb.compunds.split(' - ')
            dataToInsert.append((pdbName, pdb.interfaceParts[0], compunds[0] if len(compunds) > 1 else compunds,
                                 pdb.interfaceParts[1], compunds[1] if len(compunds) > 1 else ''))
        else:
            dataToInsert.append((pdbName, pdb.interfaceParts[0], '', pdb.interfaceParts[1], ''))

    cursor = conn.cursor()
    cursor.executemany('''
    INSERT INTO proteinComplex (PDB,UnboundChainA,NameA,UnboundChainB,NameB)
    values (?,?,?,?,?)
    ''', dataToInsert)
    conn.commit()
    conn.close()
    print('database created!')

def seperating_cavities(pdb):

    pdbName = pdb.name
    filename = pdb.file
    parser = PDBParser(PERMISSIVE=1)
    atomic_coordinates = []
    try:
        struct = parser.get_structure(pdb.file,'%s_cavities.pdb' % pdb.name)
        for model in struct:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        x,y,z = atom.get_coord()
                        atomic_coordinates.append([x,y,z])
    except:
        pass

    return atomic_coordinates



def total_hydrogen_bonds(pdbsToAnalyze):
    global PDBS_DIR, RESULTS_DIR, PARAMETERS_DIR, PARAMETERS_DIR

    conn = DBConfig.get_connection()
    cursor = conn.cursor()

    pdbsNamesToChains = dict((p[0], p[1].split(':') if len(p) > 1 else None) for p in pdbsToAnalyze)
    unfavorable_Hydroegen_bonding_residues = []
    for pdbName, chains in pdbsNamesToChains.items():
        if pdbName[0:4] not in pdbsNamesToChains: continue
        pdb_path = os.path.join(PDBS_DIR, '%s_FH.pdb' % pdbName)
        pdb = PDBReader.readFile(pdb_path, chains)
        interfaceAtoms = getInterfaceAtoms(cursor, pdb)
        bonds = hbonds(pdb)
        bonds.HDPlusDefinition = False
        cBondList = bonds.hbonds(interfaceAtoms)
        print('Calcing Hbonds for %s' % pdb.name)
        hydrogen_residue_list = []
        hydrogen_residue_dict = {}
        for donor, acceptor, eng in cBondList:
            hydrogen_residues_Dn = [(pdbName ,donor.chain, str(donor.resId), donor.residue, donor.symbol), eng]
            hydrogen_residue_list.append(hydrogen_residues_Dn)
        for donor, acceptor,eng in cBondList:
            hydrogen_residues_Ac = [(pdbName ,acceptor.chain, str(acceptor.resId), acceptor.residue, acceptor.symbol), eng]
            hydrogen_residue_list.append(hydrogen_residues_Ac)
        for k, group in groupby(sorted(hydrogen_residue_list), key=lambda x: x[0]):
            for v in group:
                if k not in hydrogen_residue_dict:
                    hydrogen_residue_dict[k] = v[1:]
                else:
                    hydrogen_residue_dict[k][0] += v[1]

        
        for key, value in hydrogen_residue_dict.items():
            if value < [-2]:
                unfavorable_Hydroegen_bonding_residues.append(key)

    return unfavorable_Hydroegen_bonding_residues                        

            
def calcEnergyTerms(pdbsToAnalyze):
    global PDBS_DIR, RESULTS_DIR, PARAMETERS_DIR, PARAMETERS_DIR
    
    output_file_list = ['Ndrieding.csv', 'interfaceVDW.csv', 'electrostatic.csv', 'electrostatic-hydrophobic.csv', 'cavity_vol.csv', 'cavity_res.csv', 'residue_depth.csv']
    if all(os.path.exists(os.path.join(RESULTS_DIR, PARAMETERS_DIR, resFile)) for resFile in output_file_list):
        print('Data already exists in result directory for energy terms.')
        return
   

    conn = DBConfig.get_connection()
    cursor = conn.cursor()
    pdbsNamesToChains = dict((p[0], p[1].split(':') if len(p) > 1 else None) for p in pdbsToAnalyze)

    for pdb, chains in pdbsNamesToChains.items():
        chains2 = "".join(chains)
        chains3=list(chains2)
        print('\nEstimating cavity volume for %s' % pdb)
        pdb_path = os.path.join(PDBS_DIR, '%s_FH.pdb' % pdb)
        pdb = PDBReader.readFile(pdb_path, chains)
        interfaceAtoms = getInterfaceAtoms(cursor, pdb)
        parser = PDBParser(PERMISSIVE=1)
        structure = parser.get_structure(pdb.name, pdb.file)
        pdb_chain_file = '%s_CR.pdb' % pdb.name
        pdb_chain_file_2 = '%s_C.pdb' % pdb.name                          
        x = PDBIO()               
        x.set_structure(structure)
        x.save('{}'.format(pdb_chain_file), ChainSelect(chains3))
        pdb = PDBParser().get_structure(pdb.name, "%s_CR.pdb" % pdb.name)
        io = PDBIO()
        io.set_structure(pdb)
        io.save('{}'.format(pdb_chain_file_2), NonHetSelect())

            

    if os.path.exists(os.path.join(RESULTS_DIR, PARAMETERS_DIR, 'cavity_volume.csv')):
        print('Skipping cavity volume calculation')
    if not os.path.exists(os.path.join(RESULTS_DIR, PARAMETERS_DIR, 'cavity_volume.csv')):
        with open(os.path.join(RESULTS_DIR,  PARAMETERS_DIR, 'cavity_volume.csv'), 'w') as cavity_file:
            print('PDB Name,Volume', file=cavity_file)
            for pdb, chains in pdbsNamesToChains.items():
                print('\nEstimating cavity volume for %s' % pdb)
                pdb_path = os.path.join('%s_C.pdb' % pdb)
                pdb = PDBReader.readFile(pdb_path, chains)
                interfaceAtoms = getInterfaceAtoms(cursor, pdb)
                cavity_volume = cavities.calculateVolume(pdb, interfaceAtoms)
                print((pdb.name + ',' + str(cavity_volume)), file=cavity_file)


    if os.path.exists(os.path.join(RESULTS_DIR, COLD_SPOTS_DIR, 'Cold_spots_residues_due_to_CH_interactions.csv')):
        print("Skipping Cold spots residues.csv")

    if not os.path.exists(os.path.join(RESULTS_DIR, COLD_SPOTS_DIR, 'Cold_spots_residues_due_to_CH_interactions.csv')):
        print("Calculating Cold spots residues.csv")
        
        unfavroable_interactions_dict = get_unfavorable_interactions(pdbsNamesToChains, cursor)

        unfavorable_residues_and_chains = []

        for pdbName, chains in pdbsNamesToChains.items():
            if pdbName[0:4] not in pdbsNamesToChains: continue
            unfavorable_interactions= unfavroable_interactions_dict[pdbName[:4]]                        
            for charged_atom, hydrophobic_atom in unfavorable_interactions:
                unfav = (pdbName, charged_atom.chain,  str(charged_atom.resId), charged_atom.residue, charged_atom.symbol)
                unfavorable_residues_and_chains.append(unfav)


        all_charged_residues_and_chains = charge_positive_negative_subbstract(pdbsNamesToChains, cursor)

        print(all_charged_residues_and_chains)
        total_pi_charge =get_pi_dicts_subtract(pdbsNamesToChains, cursor)
        print(total_pi_charge)
        total_anion_aromatic_charge = get_anion_aromatic_subtract(pdbsNamesToChains, cursor)
        Total_sub_for_cold = all_charged_residues_and_chains + total_pi_charge + total_anion_aromatic_charge
       


        final_unfavorable_residues_and_chains =  (set(unfavorable_residues_and_chains)- set(Total_sub_for_cold))


    with open(os.path.join(RESULTS_DIR, COLD_SPOTS_DIR, 'Cold_spots_residues_due_to_CH_interactions.csv'),mode = 'w') as cold_file:
        
        
        total_hydrohen_bond_energy = total_hydrogen_bonds(pdbsToAnalyze)


        final_cold_spots = []
        f_cold_spot = list((set(final_unfavorable_residues_and_chains) - set(total_hydrohen_bond_energy)))
        for one in f_cold_spot:
            final_cold_spots.append(one[:3])
        final_only_cold_spots = [', '.join(map(str, x)) for x in set(final_cold_spots)]
        print('No of cold spots due to Charged Hydrophobic interactions :', len(final_only_cold_spots), file=cold_file)
        print('PDB,Chain,Residue_number', file=cold_file)
        for two in final_only_cold_spots:
            print(two, file=cold_file)


      
    for pdbName, chains  in pdbsNamesToChains.items():
        if pdbName[0:4] not in pdbsNamesToChains: continue
        pdb_path = '%s_cavities.pdb' % pdbName[0:4]
        pdb = PDBReader.readFile(pdb_path)
        diffrent_files =seperating_cavities(pdb)
        if diffrent_files:  
            db = DBSCAN(eps=2, min_samples=10).fit(diffrent_files)
            with open(os.path.join(RESULTS_DIR, COLD_SPOTS_DIR, 'No_of_cold_spots_due_to_cavities.csv'), mode= 'w') as cluster_file:
                
                labels = db.labels_
                clusters = Counter(labels)
                total_volume = sum(Counter(labels).values())
                total_cavities = len(Counter(labels))
                for x, y in clusters.items():
                    if x == -1 :
                        total_cavities_2 = total_cavities - 1
                        print('No of Cold spots due to cavities :', total_cavities_2, file=cluster_file)
                    elif x == 0:
                        print('No of Cold spots due to cavities :', total_cavities, file=cluster_file)
                        
           

            a = np.array_split(diffrent_files,1)
            for item in a:
                b = item
            clusters = defaultdict(list)

            for i,c in enumerate(db.labels_):
                lab = b[i]

                
                I = np.concatenate([['HETATM{:^5}'.format(int(str(i)[:3])),'  O   HOH {}{:^3}  '.format(c,int(str(i)[:3])),'   {:^7.3f} {:^7.3f} {:^7.3f}'.format(lab[0],lab[1],lab[2])],[' 1.00 43.38          O']], axis=0)
                clusters[c].append(I)


            with open(os.path.join('{}Sepereated_clusters.pdb'.format(pdbName[0:4])), 'w') as cluster_file:
                for k, v in clusters.items():
                    if k != -1:
                        for x in v:
                            y = ("".join(x.tolist()))
                            cluster_file.write(y + '\n')

    
            with open(os.path.join(RESULTS_DIR, COLD_SPOTS_DIR, 'Cavities_with_the_complex.pdb'), mode='w') as open_file:
                parser = PDBParser( PERMISSIVE=1)
                pdb_io = PDBIO()
                structure1 = None
                try:
                    structure1 = parser.get_structure(pdb.file, '{}Sepereated_clusters.pdb'.format(pdbName[0:4]))
                except ValueError as e:
                    print("File is empty, There are no cold spots due to cavitites")
                structure2 = parser.get_structure(pdb.file, '%s_C.pdb' % pdbName[0:4])
                structures = [structure1, structure2]
                for struct in structures:
                    if struct:
                        pdb_io.set_structure(struct)
                        pdb_io.save(open_file)

    for pdbName, chains  in pdbsNamesToChains.items():
        if pdbName[0:4] not in pdbsNamesToChains: continue
        pdb_path = '%s_cavities.pdb' % pdbName[0:4]
        pdb = PDBReader.readFile(pdb_path)
     #   path = os.mkdir(os.path.join('results/cold_spots', '{}'.format(pdbName[0:4])))
        parser = PDBParser( PERMISSIVE=1)
        pdb_io = PDB.PDBIO()
        try:
            structure1 = parser.get_structure(pdb.file, '%sSepereated_clusters.pdb' % pdb.name)
                #structure = parser.get_structure(pdb.file, '%s.pdb' % pdb.name)


            x = []
            for chain in structure1.get_chains():
                if chain.id != '-':
                    x.append(chain.center_of_mass())
            if not x:
                x.append(1000)
            


                            
            alpha = {}
            beta = {}
            structure = parser.get_structure('c', '%s_CR.pdb' % pdb.name)
            for model in structure:
                        for chain1 in model:
                            for residue in chain1:
                                for atom in residue:
                                    if residue.get_resname() not in ["TYR" , "TRP" , "PHE",  "LYS" , "ARG" ,"MET"]:
                                        residue_id = residue.get_full_id()
                                        if atom.get_name() == 'CA':
                                            alpha[residue.get_resname() + str(residue_id[3][1])+ str(residue_id[2])] = []
                                            alpha[residue.get_resname() + str(residue_id[3][1])+ str(residue_id[2])].append(atom.get_coord())
                                        if atom.get_name() == 'CB':
                                            beta[residue.get_resname() + str(residue_id[3][1])+ str(residue_id[2])] = []
                                            beta[residue.get_resname() + str(residue_id[3][1])+ str(residue_id[2])].append(atom.get_coord())


            d = []
            x1=np.array(x)
            alpha_dic = {}   
            #beta_dic = {}
            for m, n in enumerate(x1):
                for r, k in alpha.items():
                    distance = np.linalg.norm(n - k)
                    alpha_dic[r] = []
                    alpha_dic[r].append(distance)
                #print(str(m) + ":" + str(min(alpha_dic.items(), key=lambda b: b[1])))
                    d.append(
                        {
                            'Cavity_number': m,
                            'Residue': r,
                            'distance' : distance
                    
                        })

            df=pd.DataFrame(d)

            e = []
            x1=np.array(x)
            beta_dic = {}
            for g, h in enumerate(x1):
                for s, l in beta.items():
                    distance1 = np.linalg.norm(h - l)
                    beta_dic[s] = []
                    beta_dic[s].append(distance1)
                #print(str(g)+ ":" + (str(beta_dic.items())))
                #print(str(g) + ":" + ''.join("{}: {}".format(k, v) for k, v in beta_dic.items()))
                    e.append(
                        {
                            'Cavity_number': g,
                            'Residue': s,
                            'distance' : distance1
                        })

            df1=pd.DataFrame(e)
            new_df = pd.merge(df, df1,  how='left', on=['Cavity_number','Residue'])
            new_df=new_df.replace(np.nan,0.1)

            True_False = np.where(new_df["distance_x"] > new_df["distance_y"], True, False)
            new_df["equal"] = True_False
            result= new_df.loc[new_df['equal']==True]
            result_2=result.groupby('Cavity_number', group_keys=False).apply(lambda x: x.loc[x.distance_x.idxmin()])
            result_3 = result_2.iloc[:, [1, 2]]
            #result_3.drop(columns=['distance_x'], inplace=True)
            result_3 = result_3.assign(Chain=result_3['Residue'].str[-1:], Residue_number=result_3['Residue'].str[3:-1])
            result_3.drop(columns=['Residue'], inplace=True)
            result_3.columns = ['distance', 'Chain', 'Residue_number']
            result_3.drop(columns=['distance'], inplace=True)
            #result_3.columns = ['Cavity_number', 'Residue_number', 'Cn']
            result_3.to_csv('results/cold_spots/Cold_spots_due_to_cavities.csv', sep=str(','), header=True)
        except:
            pass


        

    if os.path.exists(os.path.join(RESULTS_DIR, COLD_SPOTS_DIR, 'Cold_spots_due_to_SC_interactions.csv')):
        print('Skipping same_charge_interactions calculation')
        

    if not os.path.exists(os.path.join(RESULTS_DIR, COLD_SPOTS_DIR, 'Cold_spots_due_to_SC_interactions.csv')):
        print('calculating same_charge_interactions')

        same_charge_interactions_dict = same_charge_interface(pdbsNamesToChains, cursor)
        all_charged_residues_and_chains = charge_positive_negative_subbstract(pdbsNamesToChains, cursor)
        total_hydrohen_bond_energy = total_hydrogen_bonds(pdbsToAnalyze)

        total_favorable_charge_ions = all_charged_residues_and_chains + total_hydrohen_bond_energy


        with open(os.path.join(RESULTS_DIR, COLD_SPOTS_DIR, 'Cold_spots_due_to_SC_interactions.csv'), mode="w") as same_int:
            
            final_same_residues_and_chains = []
            for pdbName, chains in pdbsNamesToChains.items():
                if pdbName[0:4] not in pdbsNamesToChains: continue
                same_charge_pi_interactions = same_charge_interactions_dict[pdbName[:4]]
                for first_atom,  second_atom in same_charge_pi_interactions:
                    same = (pdbName, first_atom.chain, str(first_atom.resId))
                    final_same_residues_and_chains.append(same)

            s_cold_spot = list((set(final_same_residues_and_chains) - set(total_favorable_charge_ions)))
            final_same_cold_spots =[', '.join(map(str, x)) for x in s_cold_spot]
            print('No of cold spots due to Same Charge  interactions :', len(final_same_cold_spots), file=same_int)
#            print('Identified charged residues which can act as cold spots due to SC interactions', file=same_int)
            print('PDB,Chain,Residue_number', file=same_int)
            for one in final_same_cold_spots:
                print(one, file=same_int)


            os.remove(os.path.join(pdbName[0:4] + "_cavities.pdb"))
            os.remove(os.path.join("a" + pdbName[0:4]))
            os.remove(os.path.join("b" + pdbName[0:4]))

    with open(os.path.join(RESULTS_DIR, COLD_SPOTS_DIR, 'PyMOL_Representation.pml'),mode = 'w') as cold_file4:
        print("load Cavities_with_the_complex.pdb", file=cold_file4)
        print("hide lines", file=cold_file4)
        print("show cartoon", file=cold_file4)
        for pdbName, chains in pdbsNamesToChains.items():
            if len(chains[1]) == 1:
                print("select chain " + chains[0] , file=cold_file4)
                print('color purple, (sele)' , file=cold_file4)
                print("select chain " + chains[1] , file=cold_file4)
                print('color orange, (sele)' , file=cold_file4)
            else:
                print("select chain " + chains[0] , file=cold_file4)
                print('color purple, (sele)' , file=cold_file4)
                print("select chain " + chains[1][0] , file=cold_file4)
                print('color orange, (sele)' , file=cold_file4)
                print("select chain " + chains[1][1] , file=cold_file4)
                print('color green, (sele)' , file=cold_file4)

        try:
            Ch_df  = pd.read_csv(os.path.join(RESULTS_DIR, COLD_SPOTS_DIR, "Cold_spots_residues_due_to_CH_interactions.csv"), sep=',', skiprows=1, header=None)
            Ch_df.columns = ["PDB", "chain", "residue"]
            for index, row in Ch_df.iterrows():
                print('select resi ', row['residue'], "and chain ", row['chain'], file=cold_file4)
                print('show sticks, (sele)', file=cold_file4)
                print('color cyan, (sele)', file=cold_file4)
        except:
            pass
        try:
            SC_df  = pd.read_csv(os.path.join(RESULTS_DIR, COLD_SPOTS_DIR, "Cold_spots_due_to_SC_interactions.csv"), sep=',', skiprows=1, header=None)
            SC_df.columns = ["PDB", "chain", "residue"]
            for index, row in SC_df.iterrows():
                print('select resi ', row['residue'], "and chain ", row['chain'], file=cold_file4)
                print('show sticks, (sele)', file=cold_file4)
                print('color blue, (sele)', file=cold_file4)
        except:
            pass
        try:
            Cavity_df  = pd.read_csv(os.path.join(RESULTS_DIR, COLD_SPOTS_DIR, "Cold_spots_due_to_cavities.csv"), sep=',', skiprows=1, header=None)
            Cavity_df.columns = ["PDB", "chain", "residue"]
            for index, row in Cavity_df.iterrows():
                print('select resi ', row['residue'], "and chain ", row['chain'], file=cold_file4)
                print('show sticks, (sele)', file=cold_file4)
                print('color red, (sele)', file=cold_file4)
        except:
            pass



    with open(os.path.join(RESULTS_DIR, COLD_SPOTS_DIR, 'Identified_cold_spots.csv'),mode = 'w') as cold_file3:
        with open(os.path.join(RESULTS_DIR, COLD_SPOTS_DIR, 'No_of_cold_spots_due_to_cavities.csv')) as f1, open(os.path.join(RESULTS_DIR, COLD_SPOTS_DIR, 'Cold_spots_residues_due_to_CH_interactions.csv')) as f2, open(os.path.join(RESULTS_DIR, COLD_SPOTS_DIR, 'Cold_spots_due_to_SC_interactions.csv')) as f3, open(os.path.join(RESULTS_DIR, COLD_SPOTS_DIR, 'Cold_spots_due_to_cavities.csv')) as f4:
            line_file1 = f1.read()
            line_file2 = f2.read()
            line_file3 = f3.read()
            line_file4 = f4.read()
            print(line_file2,'\n',line_file3,'\n',line_file1,'\n',line_file4,'\n', file=cold_file3)

            os.remove(os.path.join(RESULTS_DIR, COLD_SPOTS_DIR, "No_of_cold_spots_due_to_cavities.csv"))
            os.remove(os.path.join(RESULTS_DIR, COLD_SPOTS_DIR, "Cold_spots_residues_due_to_CH_interactions.csv"))
            os.remove(os.path.join(RESULTS_DIR, COLD_SPOTS_DIR, "Cold_spots_due_to_SC_interactions.csv"))


    if os.path.exists(os.path.join(RESULTS_DIR, PARAMETERS_DIR, 'atom_depths.csv')):
        print('Skipping atom depth calculation')
    if not os.path.exists(os.path.join(RESULTS_DIR, PARAMETERS_DIR, 'atom_depths.csv')):
        with open(os.path.join(RESULTS_DIR, PARAMETERS_DIR, 'atom_depths.csv'), 'w') as depth:
            print('PDB,chain,resId,residue,atom symbol,depth', file=depth)
            for pdbName, chains in pdbsNamesToChains.items():
                if pdbName[0:4] not in pdbsNamesToChains: continue
                pdb_path = os.path.join(PDBS_DIR, '%s_FH.pdb' % pdbName[0:4])
                pdb = PDBReader.readFile(pdb_path, chains)
                atomDistances = residue_depth(pdb)
                for atom, distance in atomDistances:
                    print(','.join([pdbName[0:4], atom.chain, str(atom.resId), atom.residue, atom.symbol, str(distance)]), file=depth)


    cursor.close()
    conn.close()

if __name__ == "__main__":
    try:
        parser = argparse.ArgumentParser(description="Setup/download protein database based on PDB")
        parser.add_argument("pdbList", help="A file with a list of PDB to download")
        parser.add_argument("--folder", help="Name of the folder to contain downloaded files")
        parser.add_argument("dbName", help="Name of the database to create.")
        args = parser.parse_args()
        if args.dbName:
            DBConfig.DB_NAME = args.dbName
        if args.pdbList is None:
            sys.exit("Please provide a file with list of PDBs to anaylze")
        if args.dbName is None:
            sys.exit("Please provide a data base name")
        
        WORKING_DIRECTORY = args.folder if args.folder is not None else os.path.dirname(os.path.abspath(args.pdbList))
        print('WORKING DIR: %s' % WORKING_DIRECTORY)



        PDBS_DIR = os.path.join(WORKING_DIRECTORY, 'pdbs')
        pdbReader.PDBS_DIR = PDBS_DIR
        RESULTS_DIR = os.path.join(WORKING_DIRECTORY, 'results')
        new_path = (os.path.join(RESULTS_DIR, 'parameters'))
        cold_spot_path = (os.path.join(RESULTS_DIR, 'cold_spots'))
        for dir in [PDBS_DIR, RESULTS_DIR]:
            if not os.path.exists(dir):
                os.mkdir(dir)
        if not os.path.exists(new_path):
            os.makedirs(new_path)
        if not os.path.exists(cold_spot_path):
            os.makedirs(cold_spot_path)



        pdbsToAnalyzeWithChains = [pdb.strip().split("_") for pdb in open(os.path.join(args.pdbList), 'r') if
                                   pdb[0:1] != '#']  

        
        pdbsToAnalyze = [pdb[0] for pdb in pdbsToAnalyzeWithChains]
        give_permission_msms()
        give_permission_reduce()
        download_DB(pdbsToAnalyze) 
        createInterfaceCSV(pdbsToAnalyzeWithChains) 
        print('''The script will now create DB. DB is required for extra calculations
            including VDW and hydrogen bonds
            ''')
        try:
            if args.dbName:
                DBConfig.DB_NAME = args.dbName
            DBConfig.connect_and_insert()
            inserting_the_data(pdbsToAnalyzeWithChains)
           
            fillInterfacePeriphrial(pdbsToAnalyzeWithChains)
            calcEnergyTerms(pdbsToAnalyzeWithChains)

        except KeyboardInterrupt:
            print('DB will not be created. Use ./results table to see the results')
    except Exception as Argument:
        with open(os.path.join(RESULTS_DIR, COLD_SPOTS_DIR,'error.log'),'w') as error_file:
            Any_error = str(Argument)
            
            if Any_error == "A 2-dimensional array must be passed.":
                print('Please give a valid PDB ID and Chain identifiers for that specific PDB ID.', file=error_file)
            else:
                print(Any_error, file=error_file)
