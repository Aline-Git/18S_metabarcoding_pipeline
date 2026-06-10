#!/usr/bin/env python3
# -*- coding: utf-8 -*-

''' This script take as input the output cluster file from VSEARCH clustering (e.g. --cluster_smallmem)
and the original fastq file with all the sequences that were  clustered and it makes one fasta file 
per cluster with all the sequences of its members. It also create a fasta file with all the sequences
of the centroids in clusters that meet the minimum size threshold.'''


#--------------------------------------------------------------------------------------------------
#                                        Import modules
#--------------------------------------------------------------------------------------------------

import os
import sys, getopt
from io import StringIO
from os import path
import argparse
import fnmatch
import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO import FastaIO
import math


#--------------------------------------------------------------------------------------------------
#                                   Command line examples
#--------------------------------------------------------------------------------------------------

'''
~/scripts/make_cluster_files.py \
-c <cluster_uc_file> \
-f <fastx file with all the sequences> \
-o <output_folder> \
-t <input fastx type> \
-s <min_cluster_size>
 
'''

#--------------------------------------------------------------------------------------------------
#                                   Parsing command line
#--------------------------------------------------------------------------------------------------

parser = argparse.ArgumentParser(description = 'Script that creates a genome file for bedtool coverage calculation')
parser.add_argument('-c','--cluster_file', help = 'vsearch cluster output file .uc containing the info on the clusters', required = True)
parser.add_argument('-o','--output_folder', help = 'name of the output folder', required = True)
parser.add_argument('-f','--fastx_file', help = 'fastx file containing all the clustered sequences', required = True)
parser.add_argument('-t','--fastx_type', help = 'type of input fastx file (fasta or fastq)', required = True)
parser.add_argument('-s','--min_size', help = 'minimum size for a cluster to be kept', required = True)
args = parser.parse_args()

# set arguments : files and folders
cluster_filename = args.cluster_file
output_foldername = args.output_folder
fastx_filename = args.fastx_file

# set arguments : fastx type
file_type = args.fastx_type

# set arguments : minimum cluster size (5 recomended)
min_size = int(args.min_size)


#-------------------------------------------------------------------------------------------------
#                                       initialization
#--------------------------------------------------------------------------------------------------

# a dictionary to store the cluster information
d_clusters_seq = {}

readID_list_trouve = []

    
#**************************************** Program *************************************************

#--------------------------------------------------------------------------------------------------
#                1. collect the number and size information of the clusters
#--------------------------------------------------------------------------------------------------

# open the .uc cluster file and read it line by line
cluster_file = open(cluster_filename,'r')
cluster_line = cluster_file.readline()

while cluster_line : 
    
    # collect information about the clusters, lines beginning by 'C' indicates centroids related information
    split_line = cluster_line.split() 

    if split_line[0] == 'C' :
        
        # get cluster size 
        cluster_size = split_line[2]      

        # if the cluster size is larger than min_size, we collect the information
        if int(cluster_size) >= min_size :   
            
            # collect centroidID and size and store it in d_clusters_seq
            centroidID = split_line[8]
            cluster_number = split_line[1]
            d_clusters_seq[cluster_number] = [centroidID]
               
    cluster_line = cluster_file.readline()

cluster_file.close()              
        
#--------------------------------------------------------------------------------------------------
#                 3. collect the information about the reads
#--------------------------------------------------------------------------------------------------
    
        
# reopen the file to collect the information about clusters members
cluster_file = open(cluster_filename,'r')
cluster_line = cluster_file.readline()

while cluster_line : 
    # inspect only the lines begining by 'H', it indicates a cluster member
    split_line = cluster_line.split() 
    if split_line[0] == 'H' : 
        
        # get the sample ID and the centroid 
        cluster_number = split_line[1]       
        
        # we only look at the members of clusters that was already selected
        if cluster_number in d_clusters_seq : 
            
            # add the readID to the list corresponding to the cluster number
            readID = split_line[8]
            d_clusters_seq[cluster_number].append(readID)

    cluster_line = cluster_file.readline()

cluster_file.close()

# here we have a dictionnary with cluster numbers as keys and the list of sequences of the cluster as values

#--------------------------------------------------------------------------------------------------
#                      3. collect the sequences to output the fasta files
#--------------------------------------------------------------------------------------------------
    
    
print('*** writing the fasta files ***')
  
indexed_file = SeqIO.index(fastx_filename,file_type)

# for each cluster in d_clusters_seq, we will output a fasta file with all the sequences of its cluster members 
for cluster_number in d_clusters_seq : 
    fasta_file = open(output_foldername + '/ASV' + cluster_number + '.fasta','w')  
    compte = 0
    for readID in d_clusters_seq[cluster_number]:
        
        # if the sequence is in the indexed file, we add it, else a message is printed : 
        # the maximum number of sequences that we write in the output file is 1000
        if compte < 1000:
            try : 
                fasta_file.write(indexed_file[readID].format('fasta'))
                compte += 1
                readID_list_trouve.append(readID)
        
            except :
                print(readID, 'was not found')
    fasta_file.close()   
    
    # write the centroids in a dedicated centroid file
    centroid_file = open(output_foldername + '/ASV' + cluster_number + '_centroid.fasta','w')
    centroid_file.write(indexed_file[d_clusters_seq[cluster_number][0]].format('fasta'))  
    centroid_file.close()
 

