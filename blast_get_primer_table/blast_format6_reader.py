#!/usr/bin/env python
# -*- coding: utf-8 -*-

" module blast_file_reader contains the definition of the class BlastFileReader"

##########################################################################
from collections import namedtuple
##########################################################################


class blastFileReader():

    """
class BlastFileReader: can read and store information from blastn software
output generated with the option -outfmt 6 (blast+ suite)
    """

    # ------------------------------------------------------------------ #
    # Constructors/Destructors                                           #
    # ------------------------------------------------------------------ #

    def __init__(self):
        """__init__: the reader is initialized with 1 members : the dictionnary where the information will be stored """

    # Members ---------------------- #

	# self.dict[seqID] = [refID, identity%]
        self.match_list = []

        
	# a named tupple to store the current match
        self.Match = namedtuple('Match',['queryID','subjectID','identity_perc','alignment_len','query_start','query_end','subject_start','subject_end','match_sens','evalue','bitscore'])
	
    def __del__(self):
        """__del__: not implemented """
        pass

    # ------------------------------------------------------------------ #
    # Methods                                                            #
    # ------------------------------------------------------------------ #

    # public:

    def read(self, filename, minimal_length = 0, minimal_evalue = 0, key_string = None):
        """read: this function read the information in the blast_file """

        blast_file = open(filename,'r')
	
        # read the file line by line from the start
        blast_line = blast_file.readline()


        while blast_line:
	
	    # extract the information by splitting the line      
            split_line = blast_line.replace('\n','').split('\t')
            queryID = split_line[0] # extract ID of the blast query
            subjectID = split_line[1] # extract ID of the blast subject (the one in the database)
            
            identity_perc = split_line[2]
            alignment_len = split_line[3]
            
            
            query_start = split_line[6]
            query_end = split_line[7]
           
            
            subject_start = split_line[8]
            subject_end = split_line[9]
            
            # put the coordinates in order
            if int(subject_end) < int(subject_start):
                match_sens = 'reverse'
               
                temp_value = subject_end
                subject_end = subject_start
                subject_start = temp_value
            
            else :
                match_sens = 'forward'
                
                
            evalue = split_line[10]
            bitscore = split_line[11]
            
            # put these attributes in a named_tupple
            current_match = self.Match(queryID,subjectID,identity_perc,alignment_len,query_start,query_end,subject_start,subject_end,match_sens,evalue,bitscore)
            self.match_list.append(current_match)
            #print(current_match)
            
            #todo see what it does if no match
            blast_line = blast_file.readline()

        blast_file.close()  
        
        # return the list of matches satisfying the thresholds
        return self.match_list

