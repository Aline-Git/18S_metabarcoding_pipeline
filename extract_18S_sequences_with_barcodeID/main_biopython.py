#!/usr/bin/env python3
# -*- coding: utf-8 -*-

''' Given a fasta file and the output from barrnap with -k euk, bac, arc and mito, extract only the 
18S sequences and write them in a fasta file '''

###################################################################################################
#                    Extraction of the eukaryotic 18S rRNA sequence
###################################################################################################

#--------------------------------------------------------------------------------------------------
#                                        Import modules
#--------------------------------------------------------------------------------------------------

from barrnap_output_reader import barrnapOutputReader
from io import StringIO
import os
from os import path
import sys, getopt
import argparse
import fnmatch
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#--------------------------------------------------------------------------------------------------
#                                   Command line examples
#--------------------------------------------------------------------------------------------------
'''
~/scripts/extraction_of_16S-23S_IGR_from_assemblies/main.py \
  -f <input_fasta_filename> \
  -o <output_filename> \
  -b <barrnap_folder> \
  -t <fastx_type>
'''

#--------------------------------------------------------------------------------------------------
#                                   function definition
#--------------------------------------------------------------------------------------------------
       
def coordinate_overlap(match1,match2) :
    overlap = False
    
    #if match1 starts before match 2
    if min(match1.start,match1.stop) < min(match2.start,match2.stop) : 
        # there is an overlap if the minimum of match 2 is lower than the maximum of match1
        if min(match2.start,match2.stop) < max(match1.start,match1.stop) :
            overlap = True
            
    # if the first coodinate is the same, it will be overlapping anyway
    elif min(match1.start,match1.stop) == min(match2.start,match2.stop) : 
        overlap = True
        
    # if match2 starts before match1
    else :
        if min(match1.start,match1.stop) < max(match2.start,match2.stop) :
            overlap = True          

    return overlap
    
    
#--------------------------------------------------------------------------------------------------
#                                   Parsing command line
#--------------------------------------------------------------------------------------------------

parser = argparse.ArgumentParser(description = 'Script extracting the 18S sequences from fasta file given the barrnap outputs')
parser.add_argument('-f','--fastx_filename', help = 'folder where the assemblies are (fastx files)', required = True)
parser.add_argument('-b','--barrnap_output', help = 'folder where the barnap outputs are', required = True)
parser.add_argument('-o','--output_filename', help = 'name of the output folder', required = True)
parser.add_argument('-t','--type', help = 'fasta or fastq', required = True)
parser.add_argument('-s','--sampleID', help = 'sample identifier (e.g. barcodexy)',required = False)


args = parser.parse_args()

# set arguments : folders
fastx_filename = args.fastx_filename
output_filename = args.output_filename
barrnap_output = args.barrnap_output


# file type
file_type = args.type

if args.sampleID != None : 
    sampleID = '_' + args.sampleID

else : 
    sampleID = ''


#--------------------------------------------------------------------------------------------------
#                                   handling files
#--------------------------------------------------------------------------------------------------


# get the list of barnap output files
barrnap_files_list = fnmatch.filter(os.listdir(barrnap_output),'*.csv')

print(barrnap_files_list)

# test if it can be opened 
output_file = open(output_filename,'w')
output_file.close()

#**************************************** Program *************************************************

#--------------------------------------------------------------------------------------------------
#                                       initialization
#--------------------------------------------------------------------------------------------------
        

d_fastx_coord = {}

# --------------------------------------------------------------------------------------------------
#                go through the list of barrnap files make one barrnap reader per kingdom
#--------------------------------------------------------------------------------------------------

for barrnap_filename in barrnap_files_list :
    
    # extract the kingdom from the barrnap output filename ('bac', 'arc', 'euk', 'mito')
    kingdom = barrnap_filename.split('_')[-1].replace('.csv','')
    
    # initialize a reader for the euk output 
    if kingdom == 'euk':
        print('read eukaryote barrnap output')
        euk_barrnap_reader = barrnapOutputReader()
        
        # read the barrnap output file for this kingdom
        euk_matches = euk_barrnap_reader.read(barrnap_output +'/'+ barrnap_filename)
    
    # initialize a reader for the 'bact' output 
    elif kingdom == 'bact' : 
        print('read bact barrnap output')
        bact_barrnap_reader = barrnapOutputReader()
        
        # read the barrnap output file for 'bact'
        bact_matches= bact_barrnap_reader.read(barrnap_output +'/'+ barrnap_filename)
        
    # initialize a reader for the 'arc' output 
    elif kingdom == 'arc' : 
        print('read arc barrnap output')
        arc_barrnap_reader = barrnapOutputReader()
        
        # read the barrnap output file for 'arc'
        arc_matches = arc_barrnap_reader.read(barrnap_output +'/'+ barrnap_filename)
        
    # initialize a reader for the 'mito' output 
    elif kingdom == 'mito' :
        print('read mito barrnap output')
        mito_barrnap_reader = barrnapOutputReader()
        
        # read the barrnap output file for 'mito'
        mito_matches = mito_barrnap_reader.read(barrnap_output +'/'+ barrnap_filename)
    

# possible improvement : output information about the sequences other than 18S euk found   

#--------------------------------------------------------------------------------------------------
#    go through all the reads that are in at least in one dictionnary and get the coordinates and the genes
#--------------------------------------------------------------------------------------------------
   
# go through all the reads that had a match with the eukaryotic 18S rRNA sequence
for readID in euk_matches['18S_rRNA'] :
        
    match_18S = euk_matches['18S_rRNA'][readID]
    
    # check if there is a better match for other gene models
    is_match_18S = True

    # see if ther is a match with bacteria for this readID
    for gene in bact_matches :
        if readID in bact_matches[gene] :
            # check if coordinate overlap, if they do, check if the other match is better
            if coordinate_overlap(match_18S,bact_matches[gene][readID]):
                if not bact_barrnap_reader.is_better(match_18S,bact_matches[gene][readID]) :
                    is_match_18S = False
        
    # if the 18S still is the best match so far, we check the other gene models
    if is_match_18S == True : 
        
        # see if ther is a match with archaea for this readID
        for gene in arc_matches :
            if readID in arc_matches[gene] : 

                # check if coordinate overlap if they do, check if the other match is better
                if coordinate_overlap(match_18S,arc_matches[gene][readID]):
                    if not arc_barrnap_reader.is_better(match_18S,arc_matches[gene][readID]) :
                        is_match_18S = False
                            
    # check if coordinate overlap if they do, check if the other match is better                 
    if is_match_18S == True : 
        
        # see if ther is a match with archaea for this readID
        for gene in mito_matches :
            if readID in mito_matches[gene] :
                
                # check if coordinate overlap
                if coordinate_overlap(match_18S,mito_matches[gene][readID]):
                    if not mito_barrnap_reader.is_better(match_18S,mito_matches[gene][readID]) :
                        is_match_18S = False                        
                                    
         
#--------------------------------------------------------------------------------------------------
#                      extract IGR sequences and write them in a fasta file
#--------------------------------------------------------------------------------------------------

    # if the 18S is better than all the other gene models, we will extract the sequence and copy it in the output file
    if is_match_18S : 
        # construct the coordinates for samtools faidx and store them in d_fastx_coord
        d_fastx_coord[match_18S.queryID] = (match_18S.start, match_18S.stop)
         
         
# open the fastx file with Biopython
indexed_file = SeqIO.index(fastx_filename,file_type)

# create the output file for the 18S rRNA sequences detected above
output_file = open(output_filename,'w')

# we count the number of sequences selected
count=0

# we go through the sequences of the fastx file
for recordID in indexed_file :

    # if the current sequence was detected as 18S we will copy it from start coordinate to stop coordinate detected by barrnap
    if recordID in d_fastx_coord : 

        if file_type == 'fastq' :
            
            # the quality annotations are also included if the input file is a fastq file
            new_letter_annotations = {}
            new_letter_annotations['phred_quality'] = indexed_file[recordID].letter_annotations['phred_quality'][d_fastx_coord[recordID][0]-1:d_fastx_coord[recordID][1]]

            # initialize a new record for the 18S sequence
            selected_record = SeqRecord(
                indexed_file[recordID].seq[d_fastx_coord[recordID][0]-1:d_fastx_coord[recordID][1]],
                id=recordID +':' + str(d_fastx_coord[recordID][0]) + '-' + str(d_fastx_coord[recordID][1]) + sampleID,
                name=recordID +':' + str(d_fastx_coord[recordID][0]) + '-' + str(d_fastx_coord[recordID][1]) + sampleID,
                description="identified 18S_rRNA seqence selected with barrnap v0.9",
                letter_annotations = new_letter_annotations,
            )
            
            # increase the count of selected 18S sequences            
            count +=1
        
        # if the input file is a fasta we do not add a phred score
        elif file_type == 'fasta' :
            # initialize a new record for the 18S sequence
            selected_record = SeqRecord(
                indexed_file[recordID].seq[d_fastx_coord[recordID][0]-1:d_fastx_coord[recordID][1]],
                id=recordID +':' + str(d_fastx_coord[recordID][0]) + '-' + str(d_fastx_coord[recordID][1]) + sampleID,
                name=recordID +':' + str(d_fastx_coord[recordID][0]) + '-' + str(d_fastx_coord[recordID][1]) + sampleID,
                description="identified 18S_rRNA seqence selected with barrnap v0.9",
            )
            # increase the count of selected 18S sequences     
            count += 1    
            
        else : 
            print('Unknown file type. Program not working')

        output_file.write(selected_record.format(file_type))
    
    
#--------------------------------------------------------------------------------------------------
#             to improve the script : output a table with the summary of the results
#--------------------------------------------------------------------------------------------------

print(count ,' sequences were selected')

print('******** end of program  ************')












