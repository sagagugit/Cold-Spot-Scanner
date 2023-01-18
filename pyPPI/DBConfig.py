"""Configuration file for access database
"""
import sqlite3
import pkg_resources
import subprocess
import sys
from getpass import getpass
import os
import csv
from pyPPI.pdbReader import PDBReader



DB_NAME = 'x'

def get_connection():
    """
    Get connection object to local database
    """
    return sqlite3.connect(DB_NAME)


def connect_and_insert():
    """Loads teh computations to a new database
    :param pdbsToAnalyzeWithChains:
    """
    con =sqlite3.connect(DB_NAME)
    cur = con.cursor()
    con.commit()
    con.close()
