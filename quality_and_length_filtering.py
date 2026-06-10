#!/usr/bin/env python3
# -*- coding: utf-8 -*-



###################################################################################################
#                    selection of sequences above quality and length threshold given in input
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
/home/aline/scripts/quality_and_length_filtering.py \
-i <input_fastq> \
-q min_qual \
-l min_len
'''

#--------------------------------------------------------------------------------------------------
#                                   Parsing command line
#--------------------------------------------------------------------------------------------------

parser = argparse.ArgumentParser(description = 'Script that sort a fastq by quality and output a table with the seqID and their quality')
parser.add_argument('-i','--input_fastq', help = 'input fastq to sort by quality', required = True)
parser.add_argument('-l','--min_len', help = 'minimum length to keep a sequence', required = True)
parser.add_argument('-q','--min_qual', help = 'minimum quality score to keep a sequence', required = False)

args = parser.parse_args()

# set arguments : folders
input_fastq_filename = args.input_fastq

min_qual = args.min_qual
min_len = args.min_len

output_fastq_filename = input_fastq_filename.replace('.fastq','_q' + min_qual + '_l'+ min_len +'.fastq')
output_quality_filename = input_fastq_filename.replace('.fastq','_quality.csv')
output_accuracy_filename = input_fastq_filename.replace('.fastq','_accuracy.csv')

# set thresholds as floats
min_qual = float(min_qual)
min_len = float(min_len)

# Read sequences and compute average quality
records = []
for record in SeqIO.parse(input_fastq_filename, "fastq"):

    quals = np.array(record.letter_annotations["phred_quality"])
    avg_proba = np.mean(np.power(10,quals/-10))
    #add average quality and average accuracy to records
    records.append((-10*np.log10(avg_proba), record,avg_proba))

num_seq_input = len(records)


num_seq_kept = 0
# write the sequences meeting the thresholds in the output file
with open(output_fastq_filename, "w") as output_handle:
    for record in records :
        #print(type(record[0])) 
        if record[0] >= min_qual and len(record[1].seq) >= min_len  :
            num_seq_kept += 1
            SeqIO.write(record[1], output_handle, "fastq")

print(f"Sorted FASTQ file saved as {output_fastq_filename}")
print(f"{num_seq_kept} of {num_seq_input} sequences met both quality and length thresholds")

# now write the quality in a csv file

output_quality_table = open(output_quality_filename,'w')
output_accuracy_table = open(output_accuracy_filename,'w')

for record in records : 
    output_quality_table.write(record[1].id + '\t' + str(record[0]) + '\n')
    output_accuracy_table.write(record[1].id + '\t' + str(record[2]) + '\n')

output_quality_table.close()
output_accuracy_table.close()







