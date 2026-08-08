#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This script blast the primers sequences in the primer file (-f) and select extracts the results
in a summary file
"""

###################################################################################################
#   extract coordinates of construct components from blast output 6 
###################################################################################################
#--------------------------------------------------------------------------------------------------
#                                        Import modules
#--------------------------------------------------------------------------------------------------

from blast_format6_reader import blastFileReader 
from io import StringIO 
import os
from os import path 
import sys, getopt
import argparse 
import Bio
import subprocess 
from Bio import SeqIO 
from Bio.SeqRecord import SeqRecord
from Bio import SeqFeature
from Bio.SeqFeature import SeqFeature, FeatureLocation

#--------------------------------------------------------------------------------------------------
#                                   Command line examples
#--------------------------------------------------------------------------------------------------


'''
~/scripts/blast_get_primer_table/main.py \
-i INPUT_FASTA \
-f /data/VITAE/primer_set_Jamy_2019.fasta \
-p 95 \
-l 90 \
-o OUTPUT_FOLDER

'''


#--------------------------------------------------------------------------------------------------
#                                   Parsing command line
#--------------------------------------------------------------------------------------------------

parser = argparse.ArgumentParser(description = '')
parser.add_argument('-i','--input_file', help = 'input file with the raw reads in fasta format (ref)', required = True)
parser.add_argument('-f','--feature_file', help = 'file containing the primer sequences (query)', required = True)
parser.add_argument('-o','--output_folder', help = 'name of the output folder', required = True)
parser.add_argument('-p','--identity_perc', help = 'minimum identity threshold for a match to be counted', required = True)
parser.add_argument('-l','--min_perc_len', help = 'minimum percentage of the length of the sequence for a match to be counted', required = True)
parser.add_argument('-s','--skip_blast', help = 'if the blast output is already here provide the path folder to the results', required = False)

args = parser.parse_args()

# set arguments : folders
output_folder_name = args.output_folder

# set arguments : input filenames
input_filename = args.input_file
feature_filename = args.feature_file

# parameters 
identity_perc = args.identity_perc
min_perc_len = int(args.min_perc_len)
   
# set arguments : the option of skipping the blast steps : 
skip_blast = False  
if args.skip_blast != None :
    skip_blast = True
    
    # give the blast output filename, else it will be created just before the blast
    blast_output_filename = args.skip_blast


#--------------------------------------------------------------------------------------------------
#                                   handling files
#--------------------------------------------------------------------------------------------------

# build the name of the summary output table and test if it can be opened 
output_summary_name = output_folder_name + '/blastn_summary.csv'

# test the path of the output file
output_file = open(output_summary_name,'w')
output_file.close()


# folder for the output files for the blast, create if not exist
blast_folder_name = output_folder_name + '/blast'

if not os.path.exists(blast_folder_name):
    os.makedirs(blast_folder_name)     

if not os.path.exists(blast_folder_name + '/blast_db'):
    os.makedirs(blast_folder_name + '/blast_db') 
    
if not os.path.exists(blast_folder_name + '/blast_files'):
    os.makedirs(blast_folder_name + '/blast_files') 

    
#**************************************** Program *************************************************

# -------------------------------------------------------------------------------------------------
#                                       initialization
#--------------------------------------------------------------------------------------------------

#initialize a dictionary to store the feature names and length
d_feature = {}

d_count_feature = {}

d_num_feature = {}

list_num_reverse = []

#--------------------------------------------------------------------------------------------------
#       get the features names and lengths, store this information in d_feature
#--------------------------------------------------------------------------------------------------

# open the fasta file with a sequence per construct feature
records = list(SeqIO.parse(feature_filename,"fasta"))
for feature_seq in records:
    d_feature[feature_seq.id] = len(feature_seq.seq)

#--------------------------------------------------------------------------------------------------
#          blast the features on the raw reads
#--------------------------------------------------------------------------------------------------

# give a name for the indexed database that will be created by blast, print it 
blast_db_file = blast_folder_name + '/blast_db/'+ os.path.splitext(input_filename)[0].split('/')[-1]

# do this steps only if the option skip_blast was not chosen
if not skip_blast :
    # call blastn subprocess to create a blast 16S database from the input 16S refenrence
    print('\n*** blast the construct features sequences in the assembly ***')
    res1 = subprocess.call(['makeblastdb', '-in' ,input_filename, '-input_type', 'fasta', '-dbtype', 'nucl', '-title', blast_db_file, '-out', blast_db_file])
    
    print('makeblastdb terminated with attribute ', res1)

    # define the output filename
    blast_output_features_filename= blast_folder_name + '/blast_files/features_blast_output.txt'
        
    #call blastn subprocess
    res2 = subprocess.call(['blastn', '-task', 'blastn-short', '-query', feature_filename, '-db', blast_db_file ,'-out', blast_output_features_filename ,  
    '-max_target_seqs', '1000000','-perc_identity', identity_perc, '-outfmt', '6'])
    
 
    print('blastn terminated with attribute ', res2)
  
#--------------------------------------------------------------------------------------------------
#                      read the results in the blast file
#--------------------------------------------------------------------------------------------------
    
# use the blast_reader to extract the results of blastn
print('*** read the blast output ***')
blast_reader = blastFileReader()

# returns a list of Match objects (named tupples)
feature_match_list = blast_reader.read(blast_output_features_filename)  

# Match(queryID,subjectID,identity_perc,alignment_len,query_start,query_end,subject_start,subject_end,match_sens,evalue,bitscore)

# on va faire un dictionnaire d_count_feature[queryID] = feature_count

for feature_match in feature_match_list : 
    # initialize the feature count if needed
    if feature_match.subjectID not in d_count_feature : 
        d_count_feature[feature_match.subjectID] = [0,0]
        
    # if the alignment length is at least min_perc_len % of the feature length, it is counted
    if int(feature_match.alignment_len) >= min_perc_len * int(d_feature[feature_match.queryID])/100:
        if '3NDf' in feature_match.queryID or 'TAReuk454FWD1' in feature_match.queryID : 
            d_count_feature[feature_match.subjectID][0] += 1     
        elif '21R' in feature_match.queryID or '22R' in feature_match.queryID :
            d_count_feature[feature_match.subjectID][1] += 1  
            
        else : 
            print('unknown primer : ', feature_match.queryID)
    #else : 
        #print('too small, alignment rejected')    

#print(d_count_feature)

#--------------------------------------------------------------------------------------------------
#                   output the raw results in a csv file
#--------------------------------------------------------------------------------------------------


print('\n******* writing results in summary output table ',output_summary_name,' *******')

# open the output file for the basic 16S blast results per sample (barcode)
output_file = open(output_summary_name,'w')

# build the header
header ='readID\tnum_forward\tnum_reverse\n'

output_file.write(header)

# go the number of forward_features to construct each line
for readID in d_count_feature :
    output_line = readID + '\t' + str(d_count_feature[readID][0]) + '\t' + str(d_count_feature[readID][1]) + '\n'

    output_file.write(output_line)
    
# close the file
output_file.close()


print('\n***** end of program *******\n') 











