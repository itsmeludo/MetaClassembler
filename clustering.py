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


import os, sys, math, random
import subprocess, multiprocessing
#import matplotlib.pyplot as plt
import numpy as np
from itertools import product
from sklearn.cluster import DBSCAN
from distances import KL, Eucl, JSD
from sequences import signature, is_valid_sequence
from scipy.spatial.distance import pdist
from scipy.stats import scoreatpercentile
#import optics
from tristan_optics import Optics


def create_clusters(signatures,id_number,min_samples=5, epsilon=0.01):
    """From a subset of read, create clusters
    
    Parameters:
    ----------
    signatures : np.array
        matrix of read signatures
    
    Return
    ------
    positions : list
        list of integer indexes
    """

    opticss = Optics(min_samples=min_samples, epsilon=epsilon, metric=JSD).fit(signatures)
    labels = opticss.cluster(epsilon_prime=epsilon)
    
    #data=prep_optics(SetOfObjects=signatures, epsilon=0.01, MinPts=3)
    #optics_run=build_optics(SetOfObjects=data, epsilon=0.01,MinPts=3, Output_file_name="/tmp/optics"+id_number+".dat")
    centroids = []
    
    #db=ExtractDBSCAN(SetOfObjects=data, epsilon_prime=0.01)
    #signatures2 = preprocessing.scale(signatures)
    #for i in range(signatures.shape[0]-1) :
        #print( JSD(signatures[i], signatures[i+1]) )
        #1.4 unsupervised classification
    #db = DBSCAN(eps=0.01, min_samples=2, metric=JSD, algorithm='auto', 
                #leaf_size=30, p=None, random_state=None).fit(signatures)
    #labels = db.labels_
    slabels = set(labels)
    # get clusterd data point
    clusters = []
    for k in slabels :
        ind = np.where( labels == k )[0]
        if k != -1 : # clustered data points
            cluster = []
            for i in ind :
                cluster.append( i )
            clusters.append( cluster ) 
    #centroids = []
    ##print(labels)
    ## for each of the clusters, compute the centroid
    for cluster in clusters:
        sub_mat = signatures[cluster,]
        distance = pdist(sub_mat, JSD)
        sum_distance = distance.sum(axis=0)
        pos = np.argmin(sum_distance)
        centroids.append( cluster[pos] )
    #print ( centroids )
    return centroids
    #print db.components_
    
    
def random_sampling(records, sampling_size, word_size, strand):
    """Create random sample from a set
    Parameters
    ----------
    records : SeqIO
    fastq reads stored in a BioPython SeqIO object
    sampling_size : int
    number of reads to select for sampling
    
    Return
    ------
    sampling_names : list    
    read names used for the clustering
    """
    read_names = [k for k in records.keys()]
    sampling_names = []
    size = len(read_names)
    visited = np.zeros(size)
    cnt = 0
    while cnt < sampling_size:
        idx = np.random.randint(size)
        if visited[idx] == 0:
            visited[idx] = 1
            seq = records.get(read_names[idx])
            if is_valid_sequence(seq, word_size, strand):
                sampling_names.append(read_names[idx])
                cnt += 1
        if np.sum(visited==0) < sampling_size-cnt:
            # problem if not enough new positions can be visited
            raise ValueError
    #for i in range(sampling_size):
        #sampling_names.append(read_names[random.randrange(1, len(records))])
    return sampling_names
    
def run_clustering(records, list_words, word_size, strand, sampling_size, id_number): 
    """Create sample and run OPTICS clustering algorithm on it
    
    Parameters
    ----------
    records : SeqIO
        fastq reads stored in a BioPython SeqIO object
    list_words : list
        words used for signatures
    word_size : int
        size of the words used for the signatures
    strand : the strand on which to compute the signature
    sampling_size : int
        number of reads to select for sampling
    
    Return
    ------
    centroids : list
        the list of the dbscan clusters
    sampling_names : list    
        read names used for the clustering
    """
    centroids = []
    
    #sampling_size=int(len(records)/int(options.read_sampling_percent)*100)
    print("Sampling {} reads".format(str(sampling_size)))
    sampling_names = random_sampling(records, sampling_size, word_size, strand)
   
    #1.3 composition profile of sampled reads  (threaded)
    signatures=np.empty((len(sampling_names), 256), dtype=np.float)
    # TBF : just remember 
    try :
        for i, read_id in enumerate(sampling_names):
            read=records.get(read_id)
            if len(read.seq) > 0:
                signatures[i] = signature(read.seq, list_words, word_size, strand)
                #print signatures.shape
                #print(read_id+" traité")
            else:
                print(read_id+" has an empty sequence")
    except :
        print("Problem computing signature, ignore this sample")
        raise
    centroids = create_clusters(signatures,id_number)
    return centroids, sampling_names, signatures
    
    

def group_selected(selected, remove_small=10):
    """ From the different sampling steps, group together the reads (methods??)
    
    Parameters
    ----------
    remove_small : int
        small clusters are not kept 
    selected : list
        precluster, each of them containing a list of tuple(signature and names)
    
    Return
    ------
    clusters : list
        the final clusters grouped together
    """
    clusters = []
    return clusters

def cutoff_for_clusters(medoid, members, per=90, min_size_clust=20):
    """Compute the cutoff for each cluster
    
    Parameters
    ----------
    medoid : signature
        the signature of the medoids
    members : list
        the signature of the other member of the cluster
    per : int
        the percentile to use
    min_size_clust : int
        minimum size require for a cluster to computer percentile
    
    Return
    ------
    cutoff : float
        the cutoff value associated to a cluster
    """
    distances = np.array(sorted([JSD(medoid, center) for center in members]))
    if len(members)+1 > min_size_clust:
        cutoff = scoreatpercentile(distances, per)
    else:
        cutoff = np.mean(distances) # ???, standard deviation ??? TODO
    return cutoff
    
def define_centroids(records, word_size, strand, sampling_size=100, 
                     bootstrap=1, cpu=1):
    """ Create a set of centroids from partial clustering and bootstrapping
    
    Parameters
    ----------
    records ; SeqIO
        reads in fastq format 
    word_size : int
        size of the word used to compute the signature
    strand :
    sampling_size : int
        for each bootstrap, how many reads did we select for sampling
    bootstrap : int
        number of bootstrap that should be performed
    cpu : int
        number of cpu to use

    Return
    ------
    centroids_names: list
        names of the fastq reads that will act as centroids
    """
    centroids = {}
    all_sampling_names = {}
    all_sampling_signatures = {}
    list_words=[]
    #optics_bootstrap_runs=[]
    
    #creation of the nammed words
    for word in (product('ACGT' ,repeat=word_size)):
        list_words.append(''.join(word))
    
    #pool = multiprocessing.Pool(processes=cpu)
    for booti in range(bootstrap):
        centers, names, signatures = run_clustering(records, list_words, word_size, strand, sampling_size,booti)
        print(centers)
        all_sampling_names[booti] = names
        all_sampling_signatures[booti] = signatures
        #optics_bootstrap_runs[booti]= Optics(50,0.5,metric="JSD").fit(signatures)
        #centroids[booti] = optics_bootstrap_runs[booti].extractDBSCAN(0.5)
        #print optics_bootstrap_runs[booti].extractDBSCAN(0.5)
        #print(optics_bootstrap_runs[booti].cluster(0.5))
    
    ##selected = []
    ##for booti in centroids:
        ##if all_sampling_names.has_key(booti):
            ##slc = []
            ##sampling_indexes = all_sampling_names[booti]
            ##sampling_signatures = all_sampling_signatures[booti]
            ##centroids_positions = centroids[booti]
            ##for pos in centroids_positions:
                ##slc.append((sampling_signatures[pos], sampling_indexes[pos]))
            ##selected.append(slc)
     
     
    #here somewhere, merge the centroids 
     
     
     
    # read results
    centroids_names = [] #set([])
    for booti in centroids:
        if all_sampling_names.has_key(booti):
            sampling_indexes = all_sampling_names[booti]
            centroids_positions = centroids[booti]
            for pos in centroids_positions:
                centroids_names.append(sampling_indexes[pos])
                #centroids_names.add(sampling_indexes[pos])
    
    #group_centroids()
    
    #print(centroids_names)
    return centroids_names
       
def aggregate(records, clusters, clusters_index, word_size, cutoff_dist, 
              starting_cutoff=5.0):
    """ Aggregate reads in a cluster, if no cluster create a new one.
    
    Parameters
    ----------
    records ; SeqIO
        reads in fastq format 
    clusters : list
        list of signatures corresponding to a cluster
    word_size : int
        size of the word to use
    cutoff_dist : list
        for each cluster, distance cutoff below which a read can be added to the
        cluster
        
    Return
    ------
    clustered_reads : dict
        a dictionary storing for each read to which cluster he belongs
    """
    read_names=records.keys()
    clustered_reads = {}
    for i, read_id in enumerate(read_names):
        read=records.get(read_id)
        try :
            if len(read.seq) > 0:
                read_signature = signature(read.seq, list_words, word_size, strand)
        except:
            continue
        distances = np.array([JSD(read_signature, cl) for cl in clusters])
        
        # True if the value is masked, otherwise False (mean we keep it)
        mask_distances = [not dist < cutoff_dist[i] 
                          for i, dist in enumerate(distances)] 
        if not (mask_distance == True).all():
            idx = np.ma.array(distances, mask=mask_distances).argmin()
            min_dist = distances[idx]
            # the read is stored in cluster idx
            clustered_reads[read_id] = idx
    return clustered_reads
    
    
