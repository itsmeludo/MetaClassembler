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


import os, sys, math, random, re, argparse, locale
#import matplotlib.pyplot as plt
from Bio import SeqIO
from sequences import signature
from clustering import define_centroids

class CustomFormatter(argparse.RawTextHelpFormatter, 
                      argparse.ArgumentDefaultsHelpFormatter):
    pass

def main():
    env = os.environ
    env['HOME'] = './'
    #forcing locale so numpy and pylab don't mess with decimal separator digit
    loc = locale.getlocale()
    locale.setlocale(locale.LC_ALL,'en_US.utf-8')
    ###########################################################'
    useH = "%prog [-i FILE] [options]"
    parser=argparse.ArgumentParser(formatter_class=CustomFormatter)
    parser.add_argument("-i","--infile", dest="filename", 
        help="fastq input file")
    parser.add_argument("-b","--boostrap", dest="boostrap_iterations", 
        help="Number of iteration for the boostrap")
    parser.add_argument("-s","--sampling", default=1, 
        dest="read_sampling_percent", 
        help="Percentage of the reads to sample for each boostrap iteration")
    parser.add_argument("-r","--repout", dest = "outdir", 
        default=os.path.join(os.getcwd(),"Metaclassembler_output"), 
        help="Output directory")
    parser.add_argument("-c","--nb_cpu", dest = "nbcpu", default=1, type=int,
        help="Number of cpu to use")
    parser.add_argument("-m","--lgMot", type=int, default=4, 
        dest="word_size", help="word size")
    parser.add_argument("-t","--strand", default="both", dest="strand", 
        help="Strand used to compute composition: forward, reverse or both")
    parser.add_argument("-d","--distance", dest="meth_dist", default="KL", 
        help=("méthode de distance entre 2 signatures :\n"
               "* KL: Kullback-Leibler\n"
               "* Eucl : Euclidienne\n"
               "* JSD : Jensen-Shannon divergence\n"))
    #parser.add_argument("-q","--qqconque", action = "store_true", dest="afaire")
    options = parser.parse_args()
    
    # no input file? exit after a short lecture
    if(not (options.filename)):
        parser.error("Input file required")
        parser.exit(status=1)        
        
    ###########################################################'
    #MAIN
    #1. Boostrap (threaded)
        #1.1 how many reads ? Bio.SeqIO.index() ?
        #1.2 sampling
        #1.3 composition profile of sampled reads
        #1.4 unsupervised classification
        #1.5 prototype centroids definition
        
    #2. Merge/fusion of centroids accross bootstrap iterations + estimation of boostrap efficiency
    #3. Composition profile of all reads (threaded)
    #4. Sorting whole library reads in prototype centroids (threaded?)
    #5. Writing output files: one fastq per centroid.
    
    
    #1. Boostrap (threaded)
    
        #1.1 how many reads ? Bio.SeqIO.index() ?
        #1.2 sampling
    records = SeqIO.index(options.filename, "fastq")
    
    centroids = define_centroids(records, options.word_size, options.strand,
                                 sampling_size=1000, bootstrap=10, 
                                 cpu=options.nbcpu)
    
    print(centroids)
    #DEVELOPMENT ALTERNATIVE STATMENT
    
    
    #1.5 prototype centroids definition
    
    #2. Merge/fusion of centroids accross bootstrap iterations + estimation of boostrap efficiency
    #3. Composition profile of all reads (threaded)
    #4. Sorting reads in prototype centroids (threaded?)
    #5. Writing output files: one fastq per centroid.
    
    locale.setlocale(locale.LC_ALL,loc)

if __name__ == "__main__":
    main()
