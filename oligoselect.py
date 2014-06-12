#! /usr/bin/env python
import sys
import fasta
import string
import array
import tempfile
import os
import os.path
import copy
from array import array

###########################################################################     
#                                                                         #     
# C O P Y R I G H T   N O T I C E                                         #     
#                                                                         #
# Filename     : oligoselect.py                                           #
# Description  : Script to find oligos of given specifications            #
# Author(s)    : Tracy K. Teal                                            #
# Organization : California Institute of Technology                       #
# Created      : September 2002                                           #  
#                                                                         #
#  Copyright (c) 2002 by:                                                 #     
#    * California Institute of Technology                                 #     
#                                                                         #     
#    All Rights Reserved.  U.S. Government Sponsorship acknowledged.      #
#                                                                         #
#                                                                         #     
# Permission is hereby granted, free of charge, to any person             #     
# obtaining a copy of this software and associated documentation files    #     
# (the "Software"), to deal in the Software without restriction,          #     
# including without limitation the rights to use, copy, modify, merge,    #     
# publish, distribute, sublicense, and/or sell copies of the Software,    #     
# and to permit persons to whom the Software is furnished to do so,       #     
# subject to the following conditions:                                    #     
#                                                                         #     
# The above copyright notice and this permission notice shall be          #     
# included in all copies or substantial portions of the Software.         #     
#                                                                         #     
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,         #     
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF      #     
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND                   #     
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS     #     
# BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN      #     
# ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN       #     
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE        #     
# SOFTWARE.                                                               #
###########################################################################




###########################################################################  
# This script is called oligoselect.py.
# It requires fasta.py be in the same directory.
# 
# If you already know what this script does and just can't remember what order
# everything goes in at the command line, here it is.  This will be repeated
# and explained in more detail below.
# 
# The input for this script should be in the form:
# ./oligoselect.py input_file output_file oligo_length blast_database gc_min gc_max
# percent_match num_matches
# 
# This script batch processes genes/ORFS and finds oligos of a length determined
# by the user.  It can be used to find oligos for many genes/ORFs at once.  
# It checks
# to make sure they are of the specified GC content, that they are unique
# in the genome database of interest, and that they are not self-annealing.
# It looks at the region towards the 5' end, preferentially choosing oligos
# that match the criteria that are closer to the 5' end.  Right now it looks
# at the region 50 base pairs from the 5' end and 50 bp from the 3' end.  
# 
# To use this script the user must 1) be able to run it on a machine that can do
# batch blast queries.  These tools can be downloaded from NCBI.  And
# 2) the user must load the database they need to have the oligos blasted against.
# This can be done by downloading the database from NCBI or elsewhere to the
# server where this script will be run and the blast analyses done.  Then follow
# blast's 'formatdb' instructions.
# 3) Python version 2.1 or later is required.
# 
# THE COMMAND LINE INPUT
# The input for this script should be in the form:
# ./oligoselect.py input_file output_file oligo_length blast_database gc_min gc_max percent_match num_matches
# 
# input_file = the file of sequence you want to find oligos for; the file
#     format is discussed below.
# output_file = the file to write the results to
# oligo_length = how long you want your oligos to be
# blast_database = the database you want to use to blast the oligos against
#    to make sure they are unique.  This database must already be blast formatted.
# gc_min = the minimum gc_content you want your oligos to have
# gc_max = the maximum gc_content you want your oligos to have
# percent_match = percent matching allowed for uniqueness  In the blast
#    analysis, the number of base pairs that match between the oligo and the
#    second best hit is found (the first best hit is the sequence itself)
#    e.g. 35bp out of a 50bp oligo would be 70 percent_match
#    The recommended value is 20.
# num_matches = the number of matches for each gene that you want returned
#
#  INPUT FILE FORMAT
#  The input file format should be tab delimited with the following
#  column format. 1: gene or ORF name 2: ascession number
#  3: start site of the gene or ORF 4: end site 5: sequence
#  If you have data in some other format just change the columns to
#  
# 
#  OUTPUT FILE FORMAT
#  1: Gene name 2: GC content
#  3: Percent match 4: oligo sequence
# 
###########################################################################  




# *******************************************************
#   Subroutine to do blast searches to check for uniqueness of oligo
# ********************************************************

def do_blast(gene, sequence, oligo_length, blast_database):
   '''
   Create a temporary file of the oligo in fasta format so it can
   be blasted against the database
   '''
   
   temp_fasta = tempfile.mktemp()
   tf = open(temp_fasta, 'w')
   tf.write('>%s\n%s' % (gene, sequence))
   tf.close()

   # Do the blast search.  Write outout to blast_temp
   os.system('blastall -p blastn -d %s -m 8 -i %s -o blast_temp' % (blast_database, temp_fasta,))

   # Clean up - remove the temporary fasta file
   os.remove(temp_fasta)

   # Analyze the output file
   blast_input = open('blast_temp', 'r')
   lines = blast_input.readlines()
   # Take the best hit (the second row in the result file, the first
   # row is the original sequence)
   # Return percent match of 0 if no hits
   if len(lines) == 1:
      return 0
   else:
      best_hit = lines[1]
   # Get the number of matches of the best hit to the target sequence
   best_hit_cols = string.split(best_hit)
   # This relys upon the blast output being consistent.  The 4th column
   # should be the number of matches.  
   matches = float(best_hit_cols[3])
   # Find the percent matching that this is
   percent = matches/oligo_length*100
   # Return the percent match
   # print "matches", matches, "\toligo_length", oligo_length, "\tpercent", percent, "\n" 
   return percent


# *******************************************************
#   Subroutine to check to see if the oligo is self-annealing
# ********************************************************

def self_annealing(seq):
   '''
   Create a matrix with the sequence on the top and the
   reverse complement on the side.  If the bases are the
   same the matrix value is 1, otherwise it is 0.  If the
   sum of a diagonal line is more than 40% of the probe
   length, the probe is substantially self complementary
   (Li and Stormo, 2001)
   '''
   
   #print "self_annealing"


   rev_seq = fasta.reverse_complement(seq)
   end = len(seq)   # the length of the oligo sequence
   max = .4 * end   # the maximum value the sum of the diagonal can be
   anneal_flag = 0  # flag used to signal self-annealing

   # create a list of list (the matrix) of the size oligo-length x oligo-length
   tmp = end * [0]
   match_matrix = end * [0]
   for i in range(0,end):
      match_matrix[i] = copy.deepcopy(tmp) 

   # Create the matrix of 0's and 1's
   for i in range(0, end):
      for j in range(0, end):
         if seq[j] == rev_seq[i]:
            match_matrix[i][j] = 1
         else:
            match_matrix[i][j] = 0

   # Find the sum of the diagonals and determine if they're > max
   for k in range(0,end):
      diag_sum = 0
      for m in range(0, end):
         for n in range(0,end):
            # sum the diagonal
            if m == n-k:
               diag_sum = match_matrix[m][n] + diag_sum
      if diag_sum > max:
         anneal_flag = 1
   if anneal_flag == 1:
      return 1
   else:
      return 0

# ******************************************************
# * Subroutine to check self-annealing subroutine
# ******************************************************

#def test_anneal():
#   test = 'CGAT'
#   self_annealing(test)
   

# *******************************************************
# * Subroutines to count number of bases in the sequence
# *******************************************************

def a_count(sequence):
   '''
   Count the number of A's in the sequence
   '''
   
   sequence = string.upper(sequence)
   a_count = 0
   for x in range(0,len(sequence)):
      if sequence[x] == 'A':
         a_count = a_count + 1
   return float(a_count)

def t_count(sequence):
   '''
   Count the number of T's in the sequence
   '''
   
   sequence = string.upper(sequence)
   t_count = 0
   for x in range(0,len(sequence)):
      if sequence[x] == 'T':
         t_count = t_count + 1
   return float(t_count)

def g_count(sequence):
   '''
   Count the number of G's in the sequence
   '''
   
   sequence = string.upper(sequence)
   g_count = 0
   for x in range(0,len(sequence)):
      if sequence[x] == 'G':
         g_count = g_count + 1
   return float(g_count)

def c_count(sequence):
   '''
   Count the number of C's in the sequence
   '''
   
   sequence = string.upper(sequence)
   c_count = 0
   for x in range(0,len(sequence)):
      if sequence[x] == 'C':
         c_count = c_count + 1
   return float(c_count)


# *******************************************************
# * Subroutine to make sure no single base exceeds 50% of
# * the probe size
# *******************************************************

def base_fifty_percent(seq):
   '''
   Check to see if any single base exceeds 50% of the probe size
   '''
   
   seq_length = int(len(seq))
   percent50 = int(.49*seq_length)

   if (a_count(seq) > percent50) or (t_count(seq) > percent50) or (g_count(seq) > percent50) or (c_count(seq) > percent50):
      return 1 


# *******************************************************
# * Subroutine to make sure the length of any contiguous A, G, C or T
# * region is not more than 25% of the probe size
# *******************************************************

def contiguous(sequence):
   '''
   Go through every segment of the oligo that is the length
   of 25% of the oligo and determine if it is all A's, T's,
   G's or C's.  If it is, there is a contiguos region that is
   25% or more of the oligo.  Since 25% is the minimum case, this
   will find regions of 25% or more.
   '''
   
   size = int(.25*len(sequence))
   length = len(sequence) - size

   for i in range(0,length):
      if (a_count(sequence[i:i+size]) == size) or (t_count(sequence[i:i+size]) == size) or (g_count(sequence[i:i+size]) == size) or (c_count(sequence[i:i+size]) == size):
         return 1


# *******************************************************
# *******************************************************
# * Subroutines for making oligos
# *******************************************************
# *******************************************************


# *******************************************************
# * Subroutine for checking oligo criteria
# *******************************************************

def check_criteria(candidate, gc_min, gc_max, percent_match, gene, database, oligo_length):
   '''
   If the gc content is correct, the sequence is unique, no single
   base exceeds 50% of the oligo size, the length of any contiguous
   A, T, G or C region is less than 25% and it's
   not self-annealing, you've found your oligo.
   '''
   
   # Find the gc content
   gc = fasta.gc_content(candidate)

   # See if contiguous region is true or false
   contig = contiguous(candidate)

   # See if any base pair makes up > 50% of the sequence
   base_fifty = base_fifty_percent(candidate)

   # See if it's self-annealing
   self = self_annealing(candidate)

   if (gc_min < gc < gc_max) \
      and not contig \
      and not base_fifty \
      and not self:

      blast_score = do_blast(gene, candidate, oligo_length, database)
         
      if blast_score < percent_match:
         return 1, gc, blast_score
      else:
         return 0, 0, 0
   else:
      return 0, 0, 0

   
# *******************************************************
# * Subroutine for making individual oligos
# *******************************************************
   

def find_oligo(gene, seq, oligo_length, gc_min, gc_max, percent_match, num_matches, database, f):
   # ****** Find the oligos! **********
   
   # Find the length of the sequence
   seq_len = len(seq)

   # Set k = to sequence length - 100 because you don't want the match
   # to be in the first 100 base pairs from the 5' end
   # Start looking for oligos at
   k = seq_len - 100

   # Set the flag, so the loop exits if a match is found, and outputs
   # NONE if no match is ever found.
   flag = 0
   
   # Go through all sets of 'oligo_length' starting from the 5' end and
   # make sure you don't go past the start site + 50 more bases
   # k is initialized at seq-len - 50
   # So, there's 50bp at each end that aren't used

   while (k > oligo_length + 50):
      # pull out a candidate oligo
      candidate = seq[k - oligo_length:k]
      found, calc_gc, calc_blast_score = check_criteria(candidate, gc_min, gc_max, percent_match, gene, database, oligo_length)

      if found:
         # OUTPUT FILE FORMAT
         # 1: Gene name 2: GC content 3: percent match 
         # 4: oligo sequence
         print "**** Success! Found an oligo for", gene, "****\n"
         f.write('%s\t%s\t%s\t%s\n' % (gene, calc_gc, calc_blast_score, seq[k-oligo_length:k],))
         flag = flag + 1
         
               
      # If the criteria is not met, try again in the next 'oligo_length' region
      if (flag < num_matches):
         k = k - 1
      else:
         break

   # If the criteria are never met, write NONE to the output file
   if flag == 0:
      print "**** Bummer :( Did not find an oligo for", gene, " ****\n" 
      f.write('%s\t***\t***\t***\tNONE\n' % (gene,))



# *******************************************************
# * Function to make all oligos
# * This function reads in the input and writes the output
# *******************************************************

def make_oligo(input_file_name, gclist, f, oligo_length, database, gc_min, gc_max, percent_match, num_matches):
   # Make the oligos
   
   # OUTPUT
   # Output the parameters used to the top of the output file
   # If you don't care what parameters you used and just want a nice
   # output file with none of this clutter, just comment out the
   # line below.
   f.write('Input file: %s\tOligo length: %s\tDatabase: %s\tMin GC: %s\tMax GC: %s\tPercent Match: %s\n' % (input_file_name, oligo_length, database, gc_min, gc_max, percent_match))



   # INPUT
   # Read in the input file
   # Each gene/ORF should be in it's own line
   # This takes the input which is a file with columns: gene name, acession
   # number, start of the gene, end of the gene, the gene sequence

   for line in gclist.readlines():

      # Get each gene from the list
      line = string.strip(line)

      # Split it into columns
      cols = string.split(line)

   # ---------- CHANGING TO A DIFFERENT INPUT FILE FORMAT ---------
   # If your input file format is different that what has been specified, you
   # can change the columns below to reflect that

      # Take the gene name
      gene = cols[0]

      # Take the accession number
      # accession = 1
      # If you have an accession number you need, note its column here
      # and comment the above line and uncomment the below line.
      # Otherwise it's 1.
      # accession = cols[1]

      
      # Take the gene sequence
      seq = cols[4]

      # Take the start site
      # start = int(cols[1])

      find_oligo(gene, seq, oligo_length, gc_min, gc_max,
                 percent_match, num_matches, database, f)

   # ------------- END CHANGES FOR INPUT FILE FORMAT --------------







# *******************************************************
# *******************************************************
#
#   The main routine
#
# *******************************************************
# *******************************************************


# Main function - use this if called from the command line


def main():
   # Main function

   # ----------------- VARIABLES -------------------------

   input_file = sys.argv[1]
   output_file = sys.argv[2]
   oligo_length = int(sys.argv[3])
   database = sys.argv[4]
   gc_min = float(sys.argv[5])
   gc_max = float(sys.argv[6])
   percent_match = float(sys.argv[7])
   num_matches = float(sys.argv[8])
   bend = 2

   # --------------- END VARIABLES ------------------------

   # Check to make sure variables are as expected

   try:
      gclist = open(input_file, 'r')
   except:
      print "\nCannot open the input file",input_file,"\nPlease check your input file name and try again.\n"
      sys.exit(1)

   try:
      f = open(output_file, 'w')
   except:
      print "\nCannot open the output file",output_file,"\nPlease check to make sure your directory permissions allow you to write files.\n"
      sys.exit(1)

   if os.path.exists(database) == 0:
      print "\nThe database",database,"does not exist.\nPlease check the filename and try again.\n"
      sys.exit(1)

   if gc_min > 1 or gc_max > 1:
      print "\nGC_min and GC_max must be decimal values (e.g. .48 and .52)\nPlease re-enter these values and try again.\n"
      sys.exit(1)

   if percent_match < 1:
      print "\npercent_match must be a value greater than 1. \nA 20% match should be entered as 20\nPlease re-enter this value and try again.\n"
      sys.exit(1)

   if num_matches < 1:
      print "\nYou are searching for",num_matches,"matches.\nYou need to search for at least 1 match.  \nPlease re-enter the value and try again.\n"
      sys.exit(1)

   make_oligo(input_file, gclist, f, oligo_length, database, gc_min, gc_max, percent_match, num_matches)


         
if __name__ == '__main__':
   main()











