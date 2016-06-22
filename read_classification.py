#!/usr/bin/python
# -*- coding: utf-8 -*-
# <one line to give the program's name and a brief idea of what it does.>
# Copyright (C) <year>  <name of author>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
from __future__ import print_function, division

__author__ = "Ludovic Mallet, Tristan Bitard-Feildel" #je me permets
__description__ = """
Classification of long metagenomes reads in compositional groups to improve 
assemlies in low coverage conditions"
"""
__version__ = 0.001 #?
__year__ = 2014
__email__ = "l.mallet@uni-muenster.de"
__institute__ = "Insitute for Evolution and Biodiversity"
__lab__ = "Evolutionary Bioinformatics"


#import matplotlib.pyplot as plt

import numpy as np
from distances import "KL", "Eucl", "JSD"
# ignore warning comping from np.log function 
warnings.filterwarnings('ignore')


def split_cluster(cluster, matdists):
    """Split a cluster into two different clusters to reduce the variances
    """
    sub1, sub2 = None, None
    new_center1, new_center2 = None, None
    return new_center1, sub1, new_center2, sub2
    
def split_cluster_hclust(cluster, matdists):
    """Split a cluster according to a hierarchical clustering
    """
    new_clusters, new_centers = None, None
    return new_clusters, new_centers

def compute_variance(cluster):
    """Compute the variance of a cluster
    """
    var = None
    mat = None
    return var, mat

def reads_classifications(centroids, records):
    """Assign each records to a centroids
    
    """
    pass

def check_classficiation(clusters):
    """ Check the classification for each cluster, compute is variance and 
    decide to split it into two or not according to the distances
    """
    pass

    