#!/usr/bin/python
# encoding: utf-8
import numpy as np
from scipy.spatial import cKDTree as _KDTree
_sqrt = np.sqrt


class KDTree():
    """
    Wrapper for KDtree for atoms
    """

    def __init__(self, atoms):
        self._atoms = np.array(atoms)
        coord_data = np.array([atom.coord for atom in atoms])
        try:
        	self._tree = _KDTree(coord_data)
        except:
        	pass

    @staticmethod
    def construct_from_data(atoms):
        """
        Constructs KDTree from list of atoms
        :param atoms: list of atoms
        :return: KDtree
        """
        tree = KDTree(atoms)
        return tree

    def findByDistance(self, query_point, distance):
        """
        Query points around query point with certain distance
        :param query_point: point to look around
        :param distance: threshold for points selection
        :return: all points within distance of sqrt(distance)
        """
        points = self._tree.query_ball_point(query_point, _sqrt(distance))
        return self._atoms[points]

    def findNearest(self, query_point, num):
        """
        Finds nearest point to a query point
        :param query_point: query point to look around
        :param num: number of nearest neighbors
        :return: points and their distances
        """
        distances, points = self._tree.query(query_point, k=num)
        if num == 1:
            return self._atoms[points], distances
        return self._atoms[points], distances
