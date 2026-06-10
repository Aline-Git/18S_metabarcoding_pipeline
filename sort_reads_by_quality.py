#!/usr/bin/env python3
# -*- coding: utf-8 -*-

" This script sort the sequences of a fastq file by decreasing order of quality in a new fastq file"

###################################################################################################
#                               SORT READS BY AVERAGE QUALITY
###################################################################################################

#--------------------------------------------------------------------------------------------------
#                                        Import modules
#--------------------------------------------------------------------------------------------------

from Bio import SeqIO
import numpy as np
from io import StringIO
import os
from os import path
import argparse

#--------------------------------------------------------------------------------------------------
#                                   Command line examples
#--------------------------------------------------------------------------------------------------

'''
~/scripts/sort_reads_by_quality.py \
-i <input_fastq> \
-q 0
'''

#--------------------------------------------------------------------------------------------------
#                                   Parsing command line
#--------------------------------------------------------------------------------------------------

parser = argparse.ArgumentParser(description = 'Script that sort a fastq by quality and output a table with the seqID and their quality')
parser.add_argument('-i','--input_fastq', help = 'input fastq to sort by quality', required = True)
parser.add_argument('-q','--min_qual', help = 'minimum quality score', required = False)

args = parser.parse_args()

# set arguments : folders
input_fastq_filename = args.input_fastq

output_fastq_filename = input_fastq_filename.replace('.fastq','_qual_sorted.fastq')
output_quality_filename = input_fastq_filename.replace('.fastq','_quality.csv')
output_accuracy_filename = input_fastq_filename.replace('.fastq','_accuracy.csv')

min_qual = float(args.min_qual)

#**************************************** Program *************************************************

#--------------------------------------------------------------------------------------------------
#               read the sequences from the input fie and compute their average quality
#--------------------------------------------------------------------------------------------------

# Read sequences and compute average quality
records = []

# read the sequences from the fastq
for record in SeqIO.parse(input_fastq_filename, "fastq"):
    
    # for each sequence, store the phred scores in a np-array,
    quals = np.array(record.letter_annotations["phred_quality"])
    
    # Transform the quality into error probability and compute the average probability of error
    avg_proba = np.mean(np.power(10,quals/-10))
    
    # store quality and probability of error 
    records.append((-10*np.log10(avg_proba), record,avg_proba))


# Sort records by average quality (decreasing order)
records.sort(reverse=True, key=lambda x: x[0])

# Write sorted records to a new FASTQ file if their average quality meets the threshold in option
with open(output_fastq_filename, "w") as output_handle:
    for record in records :
        if record[0] >= min_qual :
            SeqIO.write(record[1], output_handle, "fastq")

print(f"Sorted FASTQ file saved as {output_fastq_filename}")

# now write the quality in a csv file

# output the average qualities and error probabilities of each read in a file
output_quality_table = open(output_quality_filename,'w')
output_accuracy_table = open(output_accuracy_filename,'w')

for record in records : 
    output_quality_table.write(record[1].id + '\t' + str(record[0]) + '\n')
    output_accuracy_table.write(record[1].id + '\t' + str(record[2]) + '\n')

output_quality_table.close()
output_accuracy_table.close()







