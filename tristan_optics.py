# -*- coding: utf-8 -*-
#!/usr/bin/python
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

__author__ = "Ludovic Mallet, Tristan Bitard-Feildel"
__description__ = """
Classification of long metagenomes reads in compositional groups to improve 
assemlies in low coverage conditions"
"""
__version__ = 0.001 #?
__year__ = 2014
__email__ = "l.mallet@uni-muenster.de"
__institute__ = "Insitute for Evolution and Biodiversity"
__lab__ = "Evolutionary Bioinformatics"

import numpy as np

################################################################################
# OPTICS
# Warning:
# C'est tres inspire d'un code source trouve en ligne et pas trop teste
################################################################################

class Optics:
    """ The main Optics class
    """
    def __init__(self, min_samples, epsilon=float('inf'), metric=None ):
        """ the init function initialise the default parameters, min_pts and eps
        """
        self.epsilon = epsilon                # maximum radius to consider
        self.min_samples = min_samples    # minimum points in cluster
        self.metric = metric
    
    def fit(self, data):
        """ read data and use the appropriate methods to perform optics 
        clustering
        """
        
        n,m = data.shape
        if (n == m and self.metric == "precomputed"):
            # its a distance matrix
            X = data
        elif (self.metric == "precomputed"):
            # Error, precomputed matrix should be symetric
            raise ValueError
        elif hasattr(self.metric, '__call__'): # check if function
            # compute pairwise distance similarity
            X = np.empty((n,n))
            for i, veci in enumerate(data):
                for j in range(i+1, len(data)):
                    vecj = data[j]
                    X[i,j] = X[j,i] = self.metric(veci, vecj)
        else:
            # Error
            raise ValueError
        self.run(X)
        return self
        
    # --------------------------------------------------------------------------
    # run the OPTICS algorithm
    # --------------------------------------------------------------------------

    def run(self, X):
        """ the True function, X is a distance matrix
        TODO : That can be changed to a dataset and a distance function
        """
        self.n = X.shape[0]       
        self.cd = np.ones( self.n  )  * float('inf')       
        self.processed = np.zeros( self.n, dtype=bool )
        self.rd = np.ones( self.n  )  * float('inf')        
        self.unprocessed = range(self.n)
        self.ordered = []
        # for each unprocessed point (p)...
        
        while self.unprocessed:
            pos = self.unprocessed[ 0 ]
            
            # mark p as processed
            # find p's neighbors
            
            self.processed[pos] = True
            self.unprocessed.remove(pos)
            self.ordered.append(pos)
            
            point_neighbors = np.where( X[pos] <= self.epsilon )[0]

            # if p has a core_distance, i.e has min_samples - 1 neighbors
            # --------------------------------------------------------------------------
            # distance from a point to its nth neighbor (n = min_cluser_size)
            # --------------------------------------------------------------------------
            
            if self.cd[pos] < float('inf') : 
                core_dist = self.cd[pos]
            elif len(point_neighbors) >= self.min_samples - 1:
                sorted_neighbors = sorted([X[pos,m] for m in point_neighbors])
                self.cd[pos] = sorted_neighbors[self.min_samples - 2]
                core_dist = self.cd[pos]
            else :
                core_dist = float('inf')

            if core_dist < float('inf') :
                
                # update reachability_distance for each unprocessed neighbor
                
                seeds = []
                #self._update(point_neighbors, pos, seeds)
                for p in [m for m in point_neighbors if not self.processed[m] ] :
                    # find new reachability distance new_rd
                    # if rd is null, keep new_rd and add n to the seed list
                    # otherwise if new_rd < old rd, update rd
                    new_rd = max( self.cd[pos], X[pos,p] )
                    #print new_rd
                    if self.rd[p] == float('inf') :
                        self.rd[p] = new_rd
                        #print new_rd
                        seeds.append(p)
                    elif new_rd < self.rd[p]:
                        self.rd[p] = new_rd
                # as long as we have unprocessed neighbors...
                while(seeds):
                    # find the neighbor n with smallest reachability distance
                    seeds.sort(key=lambda n: self.rd[n])
                    n = seeds.pop(0)
                    if not self.processed[n] :
                        # mark n as processed
                        # find n's neighbors
                    
                        self.processed[n] = True
                        self.unprocessed.remove(n)
                        self.ordered.append(n)
                        
                        n_neighbors = np.where( X[n] <= self.epsilon )[0]
                    
                        # if p has a core_distance...
                        # --------------------------------------------------------------------------
                        # distance from a point to its nth neighbor (n = min_cluser_size)
                        # --------------------------------------------------------------------------

                        if self.cd[n] < float('inf') : 
                            core_dist_n = self.cd[n]
                        elif len(n_neighbors) >= self.min_samples - 1:
                            sorted_neighbors = sorted([X[n,m] for m in n_neighbors])
                            self.cd[n] = sorted_neighbors[self.min_samples - 2]
                            core_dist_n = self.cd[n]
                        else :
                            core_dist_n = float('inf')
                            
                        if core_dist_n < float('inf') :
                            # update reachability_distance for each of n's neighbors
                            #self._update(n_neighbors, n, seeds)
                            for p in [m for m in n_neighbors if not self.processed[m] ] :
                                # find new reachability distance new_rd
                                # if rd is null, keep new_rd and add n to the seed list
                                # otherwise if new_rd < old rd, update rd
                                new_rd = max( self.cd[n], X[n,p] )
                                #print new_rd
                                if self.rd[p] == float('inf') :
                                    self.rd[p] = new_rd
                                    #print new_rd
                                    seeds.append(p)
                                elif new_rd < self.rd[p]:
                                    self.rd[p] = new_rd
        # when all points have been processed
        # return the ordered list, the reachability distance and the core distance
        return self.ordered, self.rd, self.cd
        
    def cluster(self, epsilon_prime):
        """ This function create the cluster according to a second epsilon value
        """
        clusterid = 0
        labels = np.ones( len(self.ordered) ) * -1.0        
        separators = []
        i = 0
        while i < len(self.ordered) - 1:
            j = i + 1
            obi = self.ordered[i]
            obj = self.ordered[j]
            rdi = self.rd[obi] #if self.rd[this_p] else float('inf')
            rdj = self.rd[obj] #if self.rd[next_p] else float('inf')
            
            # use an upper limit to separate the clusters
            if rdi > epsilon_prime:
                separators.append( obi )
            elif rdj > epsilon_prime:
                separators.append(  obj )
                i += 1

            i += 1
        
        if separators[-1] != len(self.ordered) :
            separators.append( len(self.ordered) ) # don't forget the last one
     
        #for start, end in separators :
        clusterid = 0 
        for i in range(len(separators) - 1):
            start = separators[i]
            end = separators[i + 1]
            if end - start > self.min_samples:
                labels[self.ordered[start:end]] = clusterid
                clusterid += 1 
        return labels

    def extractDBSCAN( self, epsilon_p ) :
        """ Don't really remmember what this one is doing, I guess its 
        extracting the clusters in  a DBSCAN fashion ??
        -> This is what is was doing in the other implementations
        """
        clusterid = 0
        labels = np.ones( len(self.ordered) ) * -1.0
        for i in range(len(self.ordered)) :
            obi = self.ordered[i]
            if self.rd[ obi ] > epsilon_p :
                if self.cd[ obi ] <= epsilon_p :
                    clusterid += 1 
                    labels[ obi ] = clusterid 
                #else :
                    #labels[ obi ] = -1 
            else :
                labels[ obi ] = clusterid 
        return  labels
