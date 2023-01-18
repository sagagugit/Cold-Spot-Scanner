#!/usr/bin/env python3
import os
import pkg_resources
import subprocess

def create_interfaceDirectory(PDBS_DIR):
    
    if PDBS_DIR is None:
        path = "./chain_A_pdbs/"
    else:
        path = os.path.join(PDBS_DIR, 'chain_A_pdbs')
        
    if not os.path.exists(path):
        os.makedirs(path)
        
    if PDBS_DIR is None:
        path = "./chain_B_pdbs/"
    else:
        path = os.path.join(PDBS_DIR, 'chain_B_pdbs')
        
    if not os.path.exists(path):
        os.makedirs(path)
        
def make_xyz_dir(xyz_destination):
    path = xyz_destination
    
    if not os.path.exists(path):
        os.makedirs(path)

def create_interface_pdb(pdb, interface, PDBS_DIR):
    create_interfaceDirectory(PDBS_DIR)
    with open(os.path.join(PDBS_DIR, pdb.name + '_FH.pdb'), 'r') as original_file:
        lines = original_file.readlines()
        
    with open(os.path.join(PDBS_DIR, 'chain_A_pdbs', pdb.name + '_chainA.pdb'), 'w') as interface_file:
        serial_numbers = []
        for line in lines:
            if line[0:6] == 'ATOM  ':
                chain, serial_number= line[21:22], int(line[6:11].strip(' '))
                if chain in pdb.interfaceParts[0]:
                    print(line, file=interface_file)
                    serial_numbers.append(serial_number)
        for line in lines:
            if line[0:6] == 'CONECT' and line[6:11] in serial_numbers:
                print(line, file=interface_file)
                
    with open(os.path.join(PDBS_DIR, 'chain_B_pdbs', pdb.name + '_chainB.pdb'), 'w') as interface_file:
        serial_numbers = []
        for line in lines:
            if line[0:6] == 'ATOM  ':
                chain, serial_number= line[21:22], int(line[6:11].strip(' '))
                if chain in pdb.interfaceParts[1]:
                    print(line, file=interface_file)
                    serial_numbers.append(serial_number)
        for line in lines:
            if line[0:6] == 'CONECT' and line[6:11] in serial_numbers:
                print(line, file=interface_file)

def calculate_cavity_volumes(PDBS_DIR):
    PDB_A_DIR = os.path.join(os.path.abspath('pdbs'), 'chain_A_pdbs')
    probe_radius = 0.1
    pdb_to_xyz = 'pyPPI/3v/xyzr/pdb_to_xyzr '
    xyz_destination = 'pdbs/chain_A_pdbs/xyz_files'
    make_xyz_dir(xyz_destination)
    
    for pdb in os.listdir(PDB_A_DIR):
        pdb_path = os.path.join('pdbs', 'chain_A_pdbs', pdb)
        xyz_file = os.path.join(xyz_destination, pdb[0:11] + '.xyzr')
        subprocess.call( (pdb_to_xyz + pdb_path + ' > ' + xyz_file), shell=True)
        subprocess.call( ('pyPPI/3v/bin/Volume.exe -i ' + xyz_file + ' -p 1.5 -g 0.5'), shell=True)
