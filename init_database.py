"""Configuration file for access database
"""
import sqlite3
import sys
from getpass import getpass
import os
import csv
import argparse

DB_NAME = ''

def createDataBase():
    con =sqlite3.connect(DB_NAME)
    cur = con.cursor()
    cur.execute('''CREATE TABLE `perAtomASA` (
                    `perAtomASA_id` INTEGER NOT NULL ,
                    `PDB` TEXT NOT NULL DEFAULT '',
                    `Chain` TEXT NOT NULL DEFAULT '',
                    `Residue` TEXT NOT NULL DEFAULT '',
                    `ResId` INTEGER NOT NULL,
                    `iCode` TEXT NOT NULL DEFAULT '',
                    `Symbol` TEXT NOT NULL DEFAULT '',
                    `Atom` TEXT NOT NULL DEFAULT '',
                    `ASA` float NOT NULL,
                    `Bfactor` float NOT NULL,
                    `Seperated` char(1) NOT NULL DEFAULT '0',
                    PRIMARY KEY (`perAtomASA_id`))''')

    cur.execute('''CREATE TABLE `interfaceDist` (
                    `interfaceDist_id` INTEGER NOT NULL,
                    `PDB` TEXT NOT NULL DEFAULT '',
                    `Chains` TEXT NOT NULL DEFAULT '',
                    `Chain` TEXT NOT NULL DEFAULT '',
                    `ResId` INTEGER NOT NULL,
                    `iCode` TEXT NOT NULL DEFAULT '',
                    `Symbol` TEXT NOT NULL DEFAULT '',
                    `Atom` TEXT NOT NULL DEFAULT '',
                    `MinDist` float NOT NULL,
                    PRIMARY KEY (`interfaceDist_id`))''')

    cur.execute('''CREATE TABLE `Ndrieding` (
                    `fullDrieding_id` INTEGER NOT NULL ,
                    `PDB` TEXT NOT NULL DEFAULT '',
                    `DonorChain` TEXT NOT NULL DEFAULT '',
                    `DonorResId` INTEGER NOT NULL,
                    `DonorICode` TEXT NOT NULL DEFAULT '',
                    `DonorSymbol` TEXT NOT NULL DEFAULT '',
                    `AccChain` TEXT NOT NULL DEFAULT '',
                    `AccResId` INTEGER NOT NULL,
                    `AccICode` TEXT NOT NULL DEFAULT '',
                    `AccSymbol` TEXT NOT NULL DEFAULT '',
                    `Energy` float NOT NULL,
                    PRIMARY KEY (`fullDrieding_id`))''')

    cur.execute('''CREATE TABLE `NinterfaceAtoms` (
                    `PDB` TEXT NOT NULL DEFAULT '',
                    `Chain` TEXT NOT NULL DEFAULT '',
                    `Residue` TEXT NOT NULL DEFAULT '',
                    `ResId` INTEGER NOT NULL,
                    `iCode` TEXT NOT NULL DEFAULT '',
                    `Symbol` TEXT NOT NULL DEFAULT '',
                    `atom` TEXT NOT NULL DEFAULT '',
                    `diffASA` double NOT NULL DEFAULT '0',
                    `pk` INTEGER NOT NULL,
                    PRIMARY KEY (`pk`))''')
    cur.execute('''CREATE TABLE `electrostat` (
                    `electrostat2_id` INTEGER NOT NULL ,
                    `PDB` TEXT NOT NULL DEFAULT '',
                    `electro` float NOT NULL,
                    `pp` int NOT NULL,
                    `mm` int NOT NULL,
                    `pm` int NOT NULL,
                    PRIMARY KEY (`electrostat2_id`))''')

    cur.execute('''CREATE TABLE `interfaceVDW` (
                    `interfaceVDW_id` INTEGER NOT NULL ,
                    `PDB` TEXT NOT NULL DEFAULT '',
                    `VDV` float NOT NULL,
                    `VDVx6` float NOT NULL,
                    `ClashV` float NOT NULL,
                    `ClashS` float NOT NULL,
                    PRIMARY KEY (`interfaceVDW_id`))''')

    cur.execute('''CREATE TABLE `proteinComplex` (
                    `proteinComplex_id` INTEGER NOT NULL ,
                    `PDB` TEXT NOT NULL DEFAULT '',
                    `UnboundChainA` TEXT NOT NULL DEFAULT '',
                    `NameA` TEXT NOT NULL DEFAULT '',
                    `UnboundChainB` TEXT NOT NULL DEFAULT '',
                    `NameB` TEXT NOT NULL DEFAULT '',
                    PRIMARY KEY (`proteinComplex_id`))''')

    cur.execute('''CREATE TABLE `nRMSD` (
                    `nRMSD_id` INTEGER NOT NULL ,
                    `Complex` TEXT NOT NULL DEFAULT '',
                    `Unbound` TEXT NOT NULL DEFAULT '',
                    `Chain` TEXT NOT NULL DEFAULT '',
                    `UnboundChain` TEXT NOT NULL DEFAULT '',
                    `RMSD` float NOT NULL,
                    `iRMSD` float NOT NULL,
                    `Atoms` INTEGER NOT NULL,
                    `iAtoms` INTEGER NOT NULL,
                    PRIMARY KEY (`nRMSD_id`))''')

    cur.execute('''CREATE TABLE `interfacePeriphrial` (
                    `interfacePeriphrial_id` INTEGER NOT NULL ,
                    `PDB` TEXT NOT NULL DEFAULT '',
                    `Chain` TEXT NOT NULL DEFAULT '',
                    `ResId` INTEGER NOT NULL,
                    `Symbol` TEXT NOT NULL DEFAULT '',
                    `Peri` float NOT NULL,
                    `PropPeri` float NOT NULL,
                    PRIMARY KEY (`interfacePeriphrial_id`))''')
    cur.execute('''CREATE TABLE `donors2` (
                    `donors2_id` INTEGER NOT NULL ,
                    `Residue` TEXT NOT NULL DEFAULT '',
                    `Symbol` TEXT NOT NULL DEFAULT '',
                    `Hydrogen` TEXT NOT NULL,
                    PRIMARY KEY (`donors2_id`))''')

    cur.execute("INSERT INTO `donors2` VALUES (18,'-','N','H'),(13,'ARG','NE','HE'),(9,'ARG','NH1','1HH1'),(10,'ARG','NH1','2HH1'),(11,'ARG','NH2','1HH2'),(12,'ARG','NH2','2HH2'),(14,'ASN','ND2','1HD2'),(15,'ASN','ND2','2HD2'),(7,'CYS','SG','HG'),(17,'GLN','NE2','1HE2 '),(8,'HIS','NE2','HE2'),(4,'LYS','NZ','1HZ'),(5,'LYS','NZ','2HZ'),(6,'LYS','NZ','3HZ'),(2,'SER','OG','HG'),(3,'THR','OG1','HG1'),(16,'TRP','NE1','HE1'),(1,'TYR','OH','HH');")


    con.commit()
    con.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Setting up the database")
    parser.add_argument("dbname", help="name of the database")
    args = parser.parse_args()
    if args.dbname:
        DB_NAME = args.dbname
        createDataBase()


