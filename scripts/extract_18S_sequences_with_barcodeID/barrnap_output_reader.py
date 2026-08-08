#!/usr/bin/env python
# -*- coding: utf-8 -*-

" module containing class barrnapOutputReader. It will read barrnap ourput files (v0.9)"

##########################################################################

import os
import re
from collections import namedtuple

##########################################################################


class barrnapOutputReader():


    # ------------------------------------------------------------------ #
    # Constructors/Destructors                                           #
    # ------------------------------------------------------------------ #

    def __init__(self):

        """__init__: initiate the reader """
        

        # Members ---------------------- #
	
		
	# the information per match, collected in the file, will be store in d_values
        self.d_values = {}
        
        
	# a named tupple to store the current match
        self.Match = namedtuple('Match',['queryID','start','stop','evalue','sens', 'status','gene'])
         
    # a dictionnary with the length of all the gene models used by barrnap
        self.d_gene_len = {}
        self.d_gene_len["5S_rRNA"]=119
        self.d_gene_len["16S_rRNA"]=1585
        self.d_gene_len["23S_rRNA"]=3232
        self.d_gene_len["5_8S_rRNA"]=156
        self.d_gene_len["18S_rRNA"]=1869
        self.d_gene_len["28S_rRNA"]=2912
        self.d_gene_len["12S_rRNA"]=954
 
    def __del__(self):
        """__del__: not implemented"""
        pass

    # ------------------------------------------------------------------ #
    # Methods                                                            #
    # ------------------------------------------------------------------ #

    # public:


    def merge_matches(self, match_list):
        ''' Sometimes barrnap only detects two partial matches (one at the
        begining of the gene and the other at the end. If this is the case, 
        we would like to consider the whole gene sequence. 
        This function merges two matches on the same read, if they fullfil 
        several conditions.        
        '''
        
        # by default, we do not merge
        merge = False
                
        # In case of multiple matches on the same read and for the same gene
        if len(match_list) > 1 :
                    
            # we go through the list of matches 
            for i in range(0, len(match_list)-1):
                merge = True

                # define the current match, the next match, the gene and the readID
                current_match = match_list[i]
                next_match = match_list[i+1]
                gene = current_match.gene
                readID = current_match.queryID
                

                # check the conditions for two matches to be merged: 
                # condition 1: if we fill the gap between these two matches, it's not exceeding 
                # 105% of the gene model length. 
                if next_match.stop - current_match.start > 1.05*self.d_gene_len[gene]:
                    merge = False
                      
                else : 

                    # condition 2 : There is no match for another gene between the two parts that we want to join.
                    # we go through the list of genes detected with the model
                    for alternative_gene in self.list_small_rRNA : 
                        
                        # we consider those that differ from the gene that we wish to merge
                        if gene != alternative_gene :
                            if readID in self.d_values[alternative_gene] :
                                
                                # we go through the matches for other genes for this read
                                for alternative_matches in self.d_values[alternative_gene][readID]:

                                    # test if there is another gene inbetween the matches for the genes we want to merge
                                    if alternative_matches.start > current_match.stop and alternative_matches.stop < next_match.start : 
                                        merge = False
                        
                # if the two conditions above are fullfilled, we will merge the two matches in one
                if merge :
                      
                    # the merged gene will start at the start of the first match and will stop at the end of the second match.
                    new_start= current_match.start
                    new_stop = next_match.stop
                    
                    # we will update the evalues. The try/except is included in case of two evalues equal to zero
                    try :
                        new_evalue = current_match.evalue * next_match.evalue /(current_match.evalue + next_match.evalue)
                    except :
                        if next_match.evalue == 0 and current_match.evalue == 0 :
                            print('Merging evalues = 0')
                            new_evalue = 0
                        else:
                            print('Problem when merging evalues')

                    # store the new match caracteristics and indicate that it was merged.
                    new_match = self.Match(readID,new_start,new_stop,new_evalue,"+","merged",gene)
                    
                    # store the indices of the two matches that were merged
                    index_to_replace = [i, i+1]
                            
                    # we only merge one time per readID and gene so break here
                    break
                        
            if merge : 
                del match_list[index_to_replace[0]:index_to_replace[1]]
                match_list.append(new_match)  
            
        return match_list                 
                             

    def remove_lists(self):
        '''
        This function only keeps the first match in self.d_values for each readID
        '''

        for gene in self.d_values : 
        
            for readID in self.d_values[gene]:
                self.d_values[gene][readID] = self.d_values[gene][readID][0]

        return 1 
        
        
#--------------------------------------------------------------------------------------------------
#                                   select the best match for a reaa
#--------------------------------------------------------------------------------------------------

	
    def select_genes(self):
    
        #self.d_values[gene][readID] = [match1, match2, match3]

        # first we will aggregate all the genes that can be 
        for gene in self.d_values : 
         
            # go through the list of the readID that have a match for that gene
            for readID in self.d_values[gene]:
                
                # if there are more than a match for this gene on this read
                if len(self.d_values[gene][readID]) > 1:        
                
                    # we put the matches in two temporary lists, one for + sens and the other for - sens
                    plus_match = []
                    minus_match = []
                
                    plus_match_start = []
                    minus_match_start = []
                
                    # 1. We put the matches in the list they belong to (sign + or -) 
                    for match in self.d_values[gene][readID]:
                        if match.sens == '+' :
                            plus_match.append(match)
                            plus_match_start.append(match.start)
                        
                        elif match.sens == '-':
                            minus_match.append(match)
                            minus_match_start.append(match.start)
                        

                    # 2. Merge under conditions the matches that are in the same direction
                    plus_match = self.merge_matches(plus_match)
                    minus_match = self.merge_matches(minus_match)
                       
                
                    # 3. See which matches we want to keep in the two temporary lists
                    best_match = None
                    # get the best match in plus_match
                    for match in plus_match + minus_match: 
                    
                        # the first match will be the current best match
                        if best_match == None : 
                            best_match = match
                        
                        # use the function 'self.is_better' to decide if the current match is better than the best match
                        elif self.is_better(match, best_match):
                            best_match = match      
                 
                    # set the best match as the only match for this gene on this read
                    self.d_values[gene][readID] = [best_match] 
       
        return 1
        
        
    def is_better(self,match1,match2):

        """select_best_match between two possible matches (match1 and match2):
        takes the highest evalue unless the difference is more than 50bp, in this case take the longest"""

        result = False  # by default, match2 is better    

        # if the difference is more than 50 bp, the length is determinent
        if match1.stop - match1.start > match2.stop - match2.start + 50 :
            result = True # match1 is better than match2
            
        # if the difference is less than 50bp, 
        elif match1.stop - match1.start < match2.stop - match2.start - 50 :
            pass  # match2 is better than match1
           
        # if the length difference is less than 50bp, we determine the best match based on the evalue        
        else : 
            
            # test if match1 has a higer evalue
            if match1.evalue < match2.evalue : 
                result = True
                
            # if it's the same evalue, we select the longest
            elif match1.evalue == match2.evalue : 
                 if match1.stop - match1.start > match2.stop - match2.start :
                     result = True
        
        return result
        
        
#--------------------------------------------------------------------------------------------------
#                                   read function
#--------------------------------------------------------------------------------------------------

    def read(self,filename):

        """read: this function will read the barrnap ourput file and store the information about matches in d_values."""

        # open the file and read it line by line
        barrnap_file = open(filename,'r')
        
        barrnap_file.readline() # pass the header
        line = barrnap_file.readline()
        
        while line :

            # get the information from the current line
            readID = line.split('\t')[0]
            start =  int(line.split('\t')[3])
            stop =  int(line.split('\t')[4])
            
            strand =  line.split('\t')[6]
            evalue =  float(line.split('\t')[5])
            

            gene =  line.split('\t')[8].split(';')[0].split('=')[1]
            
            # check if the martch was tagged as partial or complete
            if 'partial' in  line.split('\t')[8]:
                status = 'partial'
            else : 
                status = 'complete'

            # store the match in a named tupple 'Match'
            #self.Match = namedtuple('Match',['queryID','start','stop','evalue','sens','status','gene'])
            current_match = self.Match(readID,start,stop,evalue,strand, status,gene)
    
            # put the information in the dictionnary d_values, first intitialize the entry if needed
            if gene not in self.d_values : 
                self.d_values[gene] = {}
                
            # if the readID is not alredy in the dictionnary for this gene, just add the  entry
            if readID not in self.d_values[gene] : 
                self.d_values[gene][readID] = [current_match]
            
            # else we need add the current match to the list of matches for this gene and this read
            else : 
                
                self.d_values[gene][readID].append(current_match)
                
                # else we do nothing 
           
            line = barrnap_file.readline()
            
        barrnap_file.close()
        
        self.list_small_rRNA = [x for x in ["5_8S_rRNA","5S_rRNA"] if x in list(self.d_values.keys()) ]
                      
        # select the best match for each gene and each readID                      
        self.select_genes() 
        
        # transform the format of the d_values dictionnary entries into 'Matches' instead of list
        self.remove_lists()
        
        # return the dictionnary with all the matches per gene and per readID
        return self.d_values


    
#--------------------------------------------------------------------------------------------------
#                                   tests (needs improvement)
#--------------------------------------------------------------------------------------------------
	
	
# test to implement
if __name__ == '__main__':
   test = read('test_barnap_filename')
 


