#!/usr/bin/env python3
# -*- coding: utf-8 -*-

''' THis script take as input the output of VSEARCH usearch_global with the option --strand both
and select the best of the two matches for each read (if there is more than one match)

'''


#--------------------------------------------------------------------------------------------------
#                                        Import modules
#--------------------------------------------------------------------------------------------------

import os
import sys, getopt
from io import StringIO
from os import path
import argparse


#--------------------------------------------------------------------------------------------------
#                                   Command line examples
#--------------------------------------------------------------------------------------------------

'''
~/scripts/select_best_strand_ub_usearch_output.py \
-i <input_file> \
-o <output_file> 
 
'''
    
    
#--------------------------------------------------------------------------------------------------
#                                   Parsing command line
#--------------------------------------------------------------------------------------------------

parser = argparse.ArgumentParser(description = 'Script that creates a genome file for bedtool coverage calculation')
parser.add_argument('-i','--input_file', help = 'name of the usearch ouput with --strand both option', required = True)
parser.add_argument('-o','--output_file', help = 'name of the new usearch output with only one match per query', required = True)

args = parser.parse_args()

# set arguments : filnames
input_filename = args.input_file
output_filename = args.output_file


#**************************************** Program *************************************************

#--------------------------------------------------------------------------------------------------
#                 read the file 2 lines by 2 lines and output the best one in the output file
#--------------------------------------------------------------------------------------------------

input_file = open(input_filename,'r')
output_file = open(output_filename,'w')

# read the first line
input_line_1 = input_file.readline()

# go through the lines of the file (two by two)
while input_line_1 : 
    # read the second line
    input_line_2 = input_file.readline()
          
    # check that they both have the same queryID
    if input_line_1.split()[0] != input_line_2.split()[0] :
        # in case the readID is different on the second line, output a message and write the first line in the output
        print('The second sequence does not match output this one and set the second as input_line_1') 
        #print(input_line_1) 
        output_file.write(input_line_1)
        
        # The second line becomes the first line
        input_line_1 = input_line_2

    # if the two lines correspond to the same readID
    else : 
        # select the one with the highest ID percentage
        if float(input_line_1.split()[2]) >= float(input_line_2.split()[2]):
            output_file.write(input_line_1)
        else : 
            output_file.write(input_line_2)
        # and read the next line to get the first line for the next read
        input_line_1 = input_file.readline()
                
input_file.close()
output_file.close()  








































