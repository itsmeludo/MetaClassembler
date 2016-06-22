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
import warnings
# ignore warning comping from np.log function 
warnings.filterwarnings('ignore')

__all__ = ["KL", "Eucl", "JSD"]

###########################################################'

def KL(mat1, mat2):
    """the  Kullback-Leibler distance
    
    Parameters
    ----------
    mat1 : np.array
        the first vector
    mat2 : np.array
        the first vector
    
    Return
    ------
    distance : float
        the distance between the two vectors
    """
    d = mat1 * np.log(mat1/mat2)
    # lorsque qu'un mot à une fréquence nulle, on génère un NaN. pour éviter 
    # la propagation des erreurs, on ne tiens pas compte de ce mot, on le met 
    # donc à 0. Ludovic Mallet 23/03/2009.
    d[np.isnan(d)] = 0  
    d[np.isinf(d)] = 0
    return (np.sum(d))
    
def Eucl(mat1, mat2):
    """the  Euclidean distance
    
    Parameters
    ----------
    mat1 : np.array
        the first vector
    mat2 : np.array
        the first vector
    
    Return
    ------
    distance : float
        the distance between the two vectors
    """
    #d = np.power(mat1-mat2,2)
    d = pow(mat1-mat2, 2)
    d[np.isnan(d)]=0#ajout le 03/03/2010
    
#Sabine Menigaud 01/03/2010, modified on 05/03/2010 Ludovic Mallet
def JSD(mat1, mat2):
    """the  Jensen-Shannon divergence
    
    Parameters
    ----------
    mat1 : np.array
        the first vector
    mat2 : np.array
        the first vector
    
    Return
    ------
    distance : float
        the distance between the two vectors
    """
    M = (mat1+mat2) / 2
    d1 = KL(mat1, M)
    d2 = KL(mat2, M)
    d = (d1/2) + (d2/2)
    #print "test calcul JSD"
    #print d
    #d = (mat1 * np.log(mat1/M))/2 + (mat2 * np.log(mat2/M))/2
    #d[np.isnan(d)]=0
    return d
