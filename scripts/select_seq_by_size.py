#!/usr/bin/env python3
# -*- coding: utf-8 -*-


'''
This script creates a new fasta/q with only the sequences longer than the threshold specified by the 
option -l. Suppress the description from the sequence name.
'''

#--------------------------------------------------------------------------------------------------
#                                        Import modules
#--------------------------------------------------------------------------------------------------
#TODO check all these modules are required

from io import StringIO
import os
from os import path
import Bio
from Bio import SeqIO
import sys, getopt
import argparse


#--------------------------------------------------------------------------------------------------
#                                   Command line examples
#--------------------------------------------------------------------------------------------------
'''
~/scripts/select_seq_by_size.py \
  -f <input_fasta_file> \
  -t fasta \
  -l 1000 \
  -o <output_fasta_file>

'''


###################################################################################################
#
#                          parse the command line
################################################################################################### 

parser = argparse.ArgumentParser(description = 'Script that creates a genome file for bedtool coverage calculation')
parser.add_argument('-o','--output', help = 'name of the present script output file', required = True)
parser.add_argument('-f','--input', help = 'name of input file', required = True)
parser.add_argument('-t','--input_format', help = 'format of input file (fasta or fastq)', required = True)
parser.add_argument('-l','--min_read_len', help = 'minimum read length to keep a sequence', required = True)
parser.add_argument('-s','--shorten_seqname', help = 'optional provide maximum seqname length', required = False)


args = parser.parse_args()

output_file_name = args.output 
input_file_name = args.input
input_format = args.input_format

min_read_len = int(args.min_read_len)

if args.shorten_seqname != None :
    max_seqname_len = int(args.shorten_seqname)
else :
    max_seqname_len = None

######################################### open files ##############################################

output_file = open(output_file_name,'w')
input_file = open(input_file_name,'r')
########################################## init ##################################################

# i will count the number of sequences in the original file, j will count the number of sequences kept
i=0
j=0

############################## build dictionnary #########################################

# parse the input file with the SeqIO module
for seq_record in SeqIO.parse(input_file, input_format):
    i+=1 
    # if the length of the current sequence is longer than the minimum length, we will write it in the output file
    if len(seq_record.seq) >= min_read_len :
	        
	# removes the sequence description in the sequence name, just comment to remove this option
        seq_record.description = ''
    
        # write the sequence in the output file
        j+=1
        output_file.write(seq_record.format(input_format))
   
# close files 
output_file.close()
input_file.close()

# print the number of sequences kept and the total initial number of sequences analyzed
print(j, ' sequences were selected out of ',i)


