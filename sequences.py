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

__all__ = ["signature"]


import re
import numpy as np
from Bio import SeqIO

#Stackoverflow
def find_overlapping(needle, haystack):
    return len(re.findall("(?=" + needle + ")", haystack))
#!Stackoverflow

def is_valid_sequence(seq, word_size, strand):
    if(strand == "reverse"):
        seq_effective=seq.reverse_complement()
    elif(strand == "both"):
        seq_effective=seq+seq.reverse_complement()
    else:
        seq=seq_effective
    #try 1: finding occurences of words
    signature= []
    nb_words=len(seq_effective)-word_size
    #print nb_words
    if nb_words>0:
        return True
    else:
        return False
        
def signature(seq, list_words, word_size, strand):
    if(strand == "reverse"):
        seq_effective=seq.reverse_complement()
    elif(strand == "both"):
        seq_effective=seq+seq.reverse_complement()
    else:
        seq=seq_effective
    #try 1: finding occurences of words
    signature= []
    nb_words=len(seq_effective)-word_size
    #print nb_words
    if nb_words>0:
        for word in (list_words):
            #print(word)
            # need to check if nb_words != 0
            sig = find_overlapping(str(word), str(seq_effective)) / nb_words
            signature.append(sig)
    else:
        raise ValueError
    return(np.array(signature))
    #try 2: reading wordlength and eating one