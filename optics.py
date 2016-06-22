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
#import math
#import json

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


#from the pre scikit https://github.com/espg/OPTICS/blob/master/OPTICS.py


# -*- coding: utf-8 -*-

###################################
##  Written by Shane Grigsby     ##
##  Email: refuge@rocktalus.com  ##
##  Date:  May 2013              ##
###################################


## Imports ##

import sys
import scipy

from sklearn.neighbors import BallTree

from sklearn.datasets.samples_generator import make_blobs
from sklearn.preprocessing import StandardScaler

## Main Class ##

class setOfObjects(BallTree):    

    """Build balltree data structure with processing index from given data in preparation for OPTICS Algorithm

    Parameters
    ----------
    data_points: array [n_samples, n_features]"""

    def __init__(self,data_points):     

        super(setOfObjects,self).__init__(data_points)

        self._n             =   len(self.data)
        self._processed     =   scipy.zeros((self._n,1),dtype=bool) ## Start all points as 'unprocessed' ##
        self._reachability  =   scipy.ones(self._n)*scipy.inf       ## Important! ##
        self._core_dist     =   scipy.ones(self._n)*scipy.nan
        self._index         =   scipy.array(range(self._n))         ## Might be faster to use a list? ##
        self._nneighbors    =   scipy.ones(self._n,dtype=int)
        self._cluster_id    =   -scipy.ones(self._n,dtype=int)      ## Start all points as noise ##
        self._is_core       =   scipy.ones(self._n,dtype=bool)
        self._ordered_list  =   []                                  ### DO NOT switch this to a hash table, ordering is important ###
        
    ## Used in prep step ##
    def _set_neighborhood(self,point,epsilon):
        self._nneighbors[point] = self.query_radius(self.data[point], epsilon, count_only=1)[0]

    ## Used in prep step ##
    def _set_core_dist(self,point,MinPts):
        self._core_dist[point]  = self.query(self.data[point],MinPts)[0][0][-1]

## Prep Method ##

### Paralizeable! ###
def prep_optics(SetofObjects,epsilon,MinPts):

    """Prep data set for main OPTICS loop

    Parameters
    ----------
    SetofObjects: Instantiated instance of 'setOfObjects' class
    epsilon: float or int
        Determines maximum object size that can be extracted. Smaller epsilons reduce run time
    MinPts: int
        The minimum number of samples in a neighborhood to be considered a core point

    Returns
    -------
    Modified setOfObjects tree structure"""
    
    for i in SetofObjects._index:
        SetofObjects._set_neighborhood(i,epsilon)
    for j in SetofObjects._index:
        if SetofObjects._nneighbors[j] >= MinPts:
            SetofObjects._set_core_dist(j,MinPts)
    print('Core distances and neighborhoods prepped for ' + str(SetofObjects._n) + ' points.')

## Main OPTICS loop ##

def build_optics(SetOfObjects,epsilon,MinPts,Output_file_name):

    """Builds OPTICS ordered list of clustering structure

    Parameters
    ----------
    SetofObjects: Instantiated and prepped instance of 'setOfObjects' class
    epsilon: float or int
        Determines maximum object size that can be extracted. Smaller epsilons reduce run time. This should be equal to epsilon in 'prep_optics'
    MinPts: int
        The minimum number of samples in a neighborhood to be considered a core point. Must be equal to MinPts used in 'prep_optics'
    Output_file_name: string
        Valid path where write access is available. Stores cluster structure""" 

    for point in SetOfObjects._index:
        if SetOfObjects._processed[point] == False:
            expandClusterOrder(SetOfObjects,point,epsilon,
                               MinPts,Output_file_name)
                               
## OPTICS helper functions; these should not be public ##

### NOT Paralizeable! The order that entries are written to the '_ordered_list' is important! ###
def expandClusterOrder(SetOfObjects,point,epsilon,MinPts,Output_file_name):
    if SetOfObjects._core_dist[point] <= epsilon:
        while not SetOfObjects._processed[point]:
            SetOfObjects._processed[point] = True
            SetOfObjects._ordered_list.append(point)
            ## Comment following two lines to not write to a text file ##
            with open(Output_file_name, 'a') as file:
                file.write((str(point) + ', ' + str(SetOfObjects._reachability[point]) + '\n'))
                ## Keep following line! ##
                point = set_reach_dist(SetOfObjects,point,epsilon)
        print('Object Found!')
    else: 
        SetOfObjects._processed[point] = True    # Probably not needed... #


### As above, NOT paralizable! Paralizing would allow items in 'unprocessed' list to switch to 'processed' ###
def set_reach_dist(SetOfObjects,point_index,epsilon):

    ###  Assumes that the query returns ordered (smallest distance first) entries     ###
    ###  This is the case for the balltree query...                                   ###
    ###  ...switching to a query structure that does not do this will break things!   ###
    ###  And break in a non-obvious way: For cases where multiple entries are tied in ###
    ###  reachablitly distance, it will cause the next point to be processed in       ###
    ###  random order, instead of the closest point. This may manefest in edge cases  ###
    ###  where different runs of OPTICS will give different ordered lists and hence   ### 
    ###  different clustering structure...removing reproducability.                   ###
    
    distances, indices = SetOfObjects.query(SetOfObjects.data[point_index],
                                            SetOfObjects._nneighbors[point_index])
    
    ## Checks to see if there more than one member in the neighborhood ##
    if scipy.iterable(distances):

        ## Masking processed values ##
        unprocessed = indices[(SetOfObjects._processed[indices] < 1)[0].T]
        rdistances = scipy.maximum(distances[(SetOfObjects._processed[indices] < 1)[0].T],SetOfObjects._core_dist[point_index])
        SetOfObjects._reachability[unprocessed] = scipy.minimum(SetOfObjects._reachability[unprocessed], rdistances)

        ### Checks to see if everything is already processed; if so, return control to main loop ##
        if unprocessed.size > 0:            
            ### Define return order based on reachability distance ###
            return sorted(zip(SetOfObjects._reachability[unprocessed],unprocessed), key=lambda reachability: reachability[0])[0][1]
        else:
            return point_index
    else: ## Not sure if this else statement is actaully needed... ##
        return point_index

## Extract DBSCAN Equivalent cluster structure ##    

# Important: Epsilon prime should be less than epsilon used in OPTICS #
def ExtractDBSCAN(SetOfObjects, epsilon_prime):      

    """Performs DBSCAN equivalent extraction for arbitrary epsilon. Can be run multiple times.

    Parameters
    ----------
    SetOfObjects: Prepped and build instance of setOfObjects
    epsilon_prime: float or int
        Must be less than or equal to what was used for prep and build steps

    Returns
    -------
    Modified setOfObjects with cluster_id and is_core attributes."""

    # Start Cluster_id at zero, incremented to '1' for first cluster 
    cluster_id = 0                           
    for entry in SetOfObjects._ordered_list:
        if SetOfObjects._reachability[entry] > epsilon_prime:
            if SetOfObjects._core_dist[entry] <= epsilon_prime:
                cluster_id += 1
                SetOfObjects._cluster_id[entry] = cluster_id
                # Two gives first member of the cluster; not meaningful, as first cluster members do not correspond to centroids #
                ## SetOfObjects._is_core[entry] = 2     ## Breaks boolean array :-( ##
            else:
                # This is only needed for compatibility for repeated scans. -1 is Noise points #
                SetOfObjects._cluster_id[entry] = -1 
        else:
            SetOfObjects._cluster_id[entry] = cluster_id
            if SetOfObjects._core_dist[entry] <= epsilon_prime:
                # One (i.e., 'True') for core points #
                SetOfObjects._is_core[entry] = 1 
            else:
                # Zero (i.e., 'False') for non-core, non-noise points #
                SetOfObjects._is_core[entry] = 0 


##### End Algorithm #####















##reproduced from https://gist.github.com/ryangomba/1724881, no licence against reuse.




#################################################################################
## POINT
#################################################################################
 
#class Point:
    
    #def __init__(self, latitude, longitude):
        
        #self.latitude = latitude
        #self.longitude = longitude
        #self.cd = None              # core distance
        #self.rd = None              # reachability distance
        #self.processed = False      # has this point been processed?
        
    ## --------------------------------------------------------------------------
    ## calculate the distance between any two points on earth
    ## --------------------------------------------------------------------------
    
    #def distance(self, point):
        
        ## convert coordinates to radians
        
        #p1_lat, p1_lon, p2_lat, p2_lon = [math.radians(c) for c in
            #self.latitude, self.longitude, point.latitude, point.longitude]
        
        #numerator = math.sqrt(
            #math.pow(math.cos(p2_lat) * math.sin(p2_lon - p1_lon), 2) +
            #math.pow(
                #math.cos(p1_lat) * math.sin(p2_lat) -
                #math.sin(p1_lat) * math.cos(p2_lat) *
                #math.cos(p2_lon - p1_lon), 2))
 
        #denominator = (
            #math.sin(p1_lat) * math.sin(p2_lat) +
            #math.cos(p1_lat) * math.cos(p2_lat) *
            #math.cos(p2_lon - p1_lon))
        
        ## convert distance from radians to meters
        ## note: earth's radius ~ 6372800 meters
        
        #return math.atan2(numerator, denominator) * 6372800
        
    ## --------------------------------------------------------------------------
    ## point as GeoJSON
    ## --------------------------------------------------------------------------
        
    #def to_geo_json_dict(self, properties=None):
        
        #return {
            #'type': 'Feature',
            #'geometry': {
                #'type': 'Point',
                #'coordinates': [
                    #self.longitude,
                    #self.latitude,
                #]
            #},
            #'properties': properties,
        #}
        
#################################################################################
## CLUSTER
#################################################################################
 
#class Cluster:
    
    #def __init__(self, points):
        
        #self.points = points
        
    ## --------------------------------------------------------------------------
    ## calculate the centroid for the cluster
    ## --------------------------------------------------------------------------
 
    #def centroid(self):
        
        #return Point(sum([p.latitude for p in self.points])/len(self.points),
            #sum([p.longitude for p in self.points])/len(self.points))
            
    ## --------------------------------------------------------------------------
    ## calculate the region (centroid, bounding radius) for the cluster
    ## --------------------------------------------------------------------------
    
    #def region(self):
        
        #centroid = self.centroid()
        #radius = reduce(lambda r, p: max(r, p.distance(centroid)), self.points)
        #return centroid, radius
        
    ## --------------------------------------------------------------------------
    ## cluster as GeoJSON
    ## --------------------------------------------------------------------------
        
    #def to_geo_json_dict(self, user_properties=None):
        
        #center, radius = self.region()
        #properties = { 'radius': radius }
        #if user_properties: properties.update(user_properties)
        
        #return {
            #'type': 'Feature',
            #'geometry': {
                #'type': 'Point',
                #'coordinates': [
                    #center.longitude,
                    #center.latitude,
                #]
            #},
            #'properties': properties,
        #}
 
#################################################################################
## OPTICS
#################################################################################
 
#class Optics:
    
    #def __init__(self, points, max_radius, min_cluster_size):
        
        #self.points = points
        #self.max_radius = max_radius                # maximum radius to consider
        #self.min_cluster_size = min_cluster_size    # minimum points in cluster
    
    ## --------------------------------------------------------------------------
    ## get ready for a clustering run
    ## --------------------------------------------------------------------------
    
    #def _setup(self):
        
        #for p in self.points:
            #p.rd = None
            #p.processed = False
        #self.unprocessed = [p for p in self.points]
        #self.ordered = []
 
    ## --------------------------------------------------------------------------
    ## distance from a point to its nth neighbor (n = min_cluser_size)
    ## --------------------------------------------------------------------------
    
    #def _core_distance(self, point, neighbors):
 
        #if point.cd is not None: return point.cd
        #if len(neighbors) >= self.min_cluster_size - 1:
            #sorted_neighbors = sorted([n.distance(point) for n in neighbors])
            #point.cd = sorted_neighbors[self.min_cluster_size - 2]
            #return point.cd
        
    ## --------------------------------------------------------------------------
    ## neighbors for a point within max_radius
    ## --------------------------------------------------------------------------
    
    #def _neighbors(self, point):
        
        #return [p for p in self.points if p is not point and
            #p.distance(point) <= self.max_radius]
            
    ## --------------------------------------------------------------------------
    ## mark a point as processed
    ## --------------------------------------------------------------------------
        
    #def _processed(self, point):
    
        #point.processed = True
        #self.unprocessed.remove(point)
        #self.ordered.append(point)
    
    ## --------------------------------------------------------------------------
    ## update seeds if a smaller reachability distance is found
    ## --------------------------------------------------------------------------
 
    #def _update(self, neighbors, point, seeds):
        
        ## for each of point's unprocessed neighbors n...
 
        #for n in [n for n in neighbors if not n.processed]:
            
            ## find new reachability distance new_rd
            ## if rd is null, keep new_rd and add n to the seed list
            ## otherwise if new_rd < old rd, update rd
            
            #new_rd = max(point.cd, point.distance(n))
            #if n.rd is None:
                #n.rd = new_rd
                #seeds.append(n)
            #elif new_rd < n.rd:
                #n.rd = new_rd
    
    ## --------------------------------------------------------------------------
    ## run the OPTICS algorithm
    ## --------------------------------------------------------------------------
 
    #def run(self):
        
        #self._setup()
        
        ## for each unprocessed point (p)...
        
        #while self.unprocessed:
            #point = self.unprocessed[0]
            
            ## mark p as processed
            ## find p's neighbors
            
            #self._processed(point)
            #point_neighbors = self._neighbors(point)
 
            ## if p has a core_distance, i.e has min_cluster_size - 1 neighbors
 
            #if self._core_distance(point, point_neighbors) is not None:
                
                ## update reachability_distance for each unprocessed neighbor
                
                #seeds = []
                #self._update(point_neighbors, point, seeds)
                
                ## as long as we have unprocessed neighbors...
                
                #while(seeds):
                    
                    ## find the neighbor n with smallest reachability distance
                    
                    #seeds.sort(key=lambda n: n.rd)
                    #n = seeds.pop(0)
                    
                    ## mark n as processed
                    ## find n's neighbors
                    
                    #self._processed(n)
                    #n_neighbors = self._neighbors(n)
                    
                    ## if p has a core_distance...
                    
                    #if self._core_distance(n, n_neighbors) is not None:
                        
                        ## update reachability_distance for each of n's neighbors
                        
                        #self._update(n_neighbors, n, seeds)
                        
        ## when all points have been processed
        ## return the ordered list
 
        #return self.ordered
        
    ## --------------------------------------------------------------------------
    
    #def cluster(self, cluster_threshold):
        
        #clusters = []
        #separators = []
 
        #i = 0
        #while i < len(self.ordered) - 1:
                        
            #this_i = i
            #next_i = i + 1
            #this_p = self.ordered[i]
            #next_p = self.ordered[next_i]
            #this_rd = this_p.rd if this_p.rd else float('infinity')
            #next_rd = next_p.rd if next_p.rd else float('infinity')
            
            ## use an upper limit to separate the clusters
            
            #if this_rd > cluster_threshold:
                #separators.append(this_i)
            #elif next_rd > cluster_threshold:
                #separators.append(next_i)
                #i += 1
            
            ## use a jump metric to separate the clusters
            ##
            ## if this_rd - next_rd > cluster_threshold:
            ##    separators.append(this_i)
            ## elif next_rd - this_rd > cluster_threshold:
            ##    separators.append(next_i)
            ##    i += 1
 
            #i += 1
                
        #for i in range(len(separators) - 1):
            #start = separators[i] + 1
            #end = separators[i + 1]
            #if end - start > self.min_cluster_size:
                #clusters.append(Cluster(self.ordered[start:end]))
 
        #return clusters
 
## LOAD SOME POINTS
 
##optics = Optics(points, 200, 5)   # 200 meter radius for neighbor consideration
##ordered = optics.run()
##clusters = optics.cluster(100)    # 100 meter threshold for clustering

