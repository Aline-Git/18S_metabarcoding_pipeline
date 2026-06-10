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
import Bio
from Bio import SeqIO
from collections import namedtuple


#--------------------------------------------------------------------------------------------------
#                                   Command line examples
#--------------------------------------------------------------------------------------------------

'''
~/scripts/make_grouped_taxonomic_table_from_usearch_output.py \
-c <cluster_uc_file> \
-b <centroid_taxonomic_table> \
-o <output_folder> 
 
'''

#--------------------------------------------------------------------------------------------------
#                                   function definition
#--------------------------------------------------------------------------------------------------
        
# this functions builds a reference dictionnary from a csv file (space or tab separated),
# with first column (id) as key and the rest as value
def build_ref_dictionnary(reference_filename):
    
    dictionary = {}
    
    # open the tax file and read it line by line
    reference_file = open(reference_filename,'r')

    #reference_file.readline() # uncomment to pass the header
    line = reference_file.readline()
    while line:
        
        # extract the referenceID and the corresponding taxonomy
        try : 
            field = line.split(' ',1)
            dictionary[field[0]]= field[1].replace('\n','').replace(';','|')    
        except : 
            field = line.split('\t',1)
            dictionary[field[0]]= field[1].replace('\n','').replace(';','|')    
        
        line = reference_file.readline() 
        
    reference_file.close()
    
    # add supplementary entries for 'No match'
    dictionary['*']='No match'
    return dictionary
    
    
#--------------------------------------------------------------------------------------------------
#                                   Parsing command line
#--------------------------------------------------------------------------------------------------

parser = argparse.ArgumentParser(description = 'Script that creates a genome file for bedtool coverage calculation')
parser.add_argument('-c','--cluster_file', help = 'vsearch cluster output file .uc containing the info on the clusters', required = True)
parser.add_argument('-o','--output_folder', help = 'name of the output folder', required = True)
parser.add_argument('-b','--tax_table', help = 'taxonomic table from the comparison of the centroids with the database', required = True)
parser.add_argument('-s','--min_size', help = 'minimum size required to consider a cluster', required = True)
parser.add_argument('-t','--corresp_file', help = 'file with correspondance between tax and accessionID', required = True)

args = parser.parse_args()

# set arguments : folders and files
cluster_filename = args.cluster_file
output_foldername = args.output_folder

# set arguments : tax_table
tax_table = args.tax_table

# csv table with the referenceID and the corresponding taxonomy
corresp_file = args.corresp_file

# set minimum cluster size
min_size = int(args.min_size)


#-------------------------------------------------------------------------------------------------
#                                       initialization
#--------------------------------------------------------------------------------------------------

# a dictionary to store the cluster information
d_cluster_info = {}

# a dictionnary to store the number of read per sample and ASV
d_ASV_table = {}

# Object Cluster with different cluster caracteristics
Cluster = namedtuple('Cluster',['cluster_number','size','taxonomy','accession_ID','identity_perc','alignment_len'])

# list of the centroids
list_centroids = []

# build the dictionary that makes a correspondance from refID to ref name
ref_to_name_dict = build_ref_dictionnary(corresp_file)
    
    
#**************************************** Program *************************************************

#--------------------------------------------------------------------------------------------------
#                1. collect the number and size information of the clusters
#--------------------------------------------------------------------------------------------------

# open the cluster file and read it line by line
cluster_file = open(cluster_filename,'r')
cluster_line = cluster_file.readline()

while cluster_line : 
    split_line = cluster_line.split() 
          
    # collect information about cluster centroids
    if split_line[0] == 'C' : 
        
        # collect centroidID and size
        cluster_size = split_line[2]
        if int(cluster_size) >= min_size :  
            centroidID = split_line[8].split(':')[0]
            cluster_size = split_line[2]
            cluster_number = split_line[1]
           
            # add an entry for the current cluster as a named tupple
            d_cluster_info[centroidID] = Cluster(cluster_number, cluster_size, '','','','')    
    
            # add a count for the centroid 
            sampleID = split_line[8].split('_',1)[1]
            
            # add an entry for the current sampleID if needed
            if sampleID not in d_ASV_table:
                d_ASV_table[sampleID] = {}
            # initialize the count with one for this cluster and this sample
            d_ASV_table[sampleID][centroidID] = 1

    cluster_line = cluster_file.readline()
cluster_file.close()

# now we reopen the file to collect the information about the clusters members
cluster_file = open(cluster_filename,'r')
cluster_line = cluster_file.readline()

while cluster_line : 
    split_line = cluster_line.split() 
    
    # collect information about the cluster members
    if split_line[0] == 'H' : 
        # get the centroidID of the cluster
        centroidID = split_line[9].split(':')[0]
        
        # we count the read only if the cluster size is at least the minimum size
        if centroidID in d_cluster_info :
            
            # get the sampleID 
            sampleID = split_line[8].split('_',1)[1]

            # add a count in the ASV table
            # initialize an entry with the current sampleID if needed
            if sampleID not in d_ASV_table : 
                d_ASV_table[sampleID] = {}
                
            # initialize an entry with the current centroid if needed
            if centroidID not in d_ASV_table[sampleID] : 
                d_ASV_table[sampleID][centroidID] = 0
            
            # add a count for this cluster member
            d_ASV_table[sampleID][centroidID] += 1
    cluster_line = cluster_file.readline()

cluster_file.close()


#--------------------------------------------------------------------------------------------------
#                      2. collect information from the taxonomic table
#--------------------------------------------------------------------------------------------------
    
# use the blast_reader to extract the results of blastn, ie assign 
print('*** reading the taxonomic table of the centroids ***')
    
# open the tax file (output from the comparison of the sequences with the database) and read it line by line
tax_file = open(tax_table,'r')
tax_line = tax_file.readline()

while tax_line : 
    split_tax = tax_line.split()
    
    # get the centroidID 
    centroidID = split_tax[0].split('_')[0]

    # get information about the match and add them to the dictionary d_cluster_info
    d_cluster_info[centroidID] = d_cluster_info[centroidID]._replace(taxonomy = ref_to_name_dict[split_tax[1]])
    d_cluster_info[centroidID] = d_cluster_info[centroidID]._replace(identity_perc = split_tax[2])
    d_cluster_info[centroidID] = d_cluster_info[centroidID]._replace(alignment_len = split_tax[3])
    d_cluster_info[centroidID] = d_cluster_info[centroidID]._replace(accession_ID = split_tax[1])
    
    tax_line = tax_file.readline()
    
  
#--------------------------------------------------------------------------------------------------
#                  3. output the cluster info table
#--------------------------------------------------------------------------------------------------
    
print('\n******* writing the ASV information in  ',output_foldername + '/ASV_info.csv',' *******')

# create a ASV_info table with information about the clusters
output_file = open(output_foldername + '/ASV_info.csv','w')

# write one output line per referenceID (key)
header = "ASV\ttaxonomy\tsize\tcentroidID\taccessionID\tIDperc\talignement_len\n"
    
output_file.write(header)

for centroidID in d_cluster_info : 
    
    # the row name will be the ASV number
    output_line = 'ASV' + d_cluster_info[centroidID].cluster_number 

    # add the percentage id and the alignment length
    output_line += "\t" + d_cluster_info[centroidID].taxonomy 
    output_line += '\t' + d_cluster_info[centroidID].size
    output_line += '\t' + centroidID
    output_line += '\t' + d_cluster_info[centroidID].accession_ID
    output_line += '\t' + d_cluster_info[centroidID].identity_perc 
    output_line += '\t' + d_cluster_info[centroidID].alignment_len + '\n'

    output_file.write(output_line)

output_file.close()


  
#--------------------------------------------------------------------------------------------------
#                  3. output the ASV table
#--------------------------------------------------------------------------------------------------
    
print('\n******* writing the ASV table in  ',output_foldername + '/ASV_table.csv',' *******')
output_file = open(output_foldername + '/ASV_table.csv','w')


# the header is ASV\tASV3\tASV420, etc
header = "ASV"
for centroidID in d_cluster_info : 
    header += '\t' + 'ASV' + d_cluster_info[centroidID].cluster_number    

header += '\n'
output_file.write(header)

for sampleID in d_ASV_table : 
#    print(sampleID , 'in output')    
    output_line = sampleID

    for centroidID in d_cluster_info : 
        
        if centroidID in d_ASV_table[sampleID] : 
            output_line += '\t' + str(d_ASV_table[sampleID][centroidID])
        else : 
            output_line += '\t' + '0'
 #       print('output_line : ',output_line)
    output_line += '\n'

    output_file.write(output_line)

output_file.close()












































