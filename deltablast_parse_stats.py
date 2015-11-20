# deltablast-stats.py
# David Shin
# 2013-11-14
# Version 
# Cut deltablast-parse6i.py to make a screener: 

import time
import os
import sys
import glob
import shutil
from Bio import SeqIO 
from Bio.Blast import NCBIXML
start = time.time()

def read_param_file():
                                   # define globals
    global high_identity, low_identity, high_length, low_length
    global input_file_names, output_file_root, output_file_path
    global w_header, w_names, w_lengths, w_identities, w_alignment
    global output_alignment_length
                                   # opens a file to read
#    from_file = open('parse6.inp', 'r')
    filename = sys.argv[-1]
    from_file = open(filename, 'r')
                                   # reads the file into a list
    my_file_list = []
    my_file_list = from_file.readlines()
                                   # assign variables & convert to float if needed
    input_file_names = my_file_list[9][:-1]

    high_identity = float(my_file_list[16]) / 100
    low_identity  = float(my_file_list[19]) / 100
    high_length   = float(my_file_list[22]) / 100  
    low_length    = float(my_file_list[25]) / 100 
 
    output_file_root = my_file_list[32][:-1]
    output_file_path = output_file_root + "_results" + "/" 

    w_header     = my_file_list[39][:-1]
    w_names      = my_file_list[42][:-1]
    w_lengths    = my_file_list[45][:-1]
    w_identities = my_file_list[48][:-1]
    w_alignment  = my_file_list[51][:-1]
    
    output_alignment_length = int(my_file_list[54][:-1]) 
                           
    from_file.close()

def check_output_dir():
    dir = output_file_path
    if not os.path.exists(dir):
        os.makedirs(dir)
    else:
        shutil.rmtree(dir)           #removes all the subdirectories!
        os.makedirs(dir)
            
def initial_alignment_length():
    global first, second 
    first = 0
    second = output_alignment_length
#    print "dali lamma" #get rid of this later
    
def write_param_file():
    high_identity2 = high_identity * 100
    low_identity2  = low_identity * 100
    high_length2   = high_length * 100
    low_length2    = low_length * 100

    fob2 = open(output_file_root + ".def", 'w')
    
    fob2.write("high percent identity cutoff: " + str(high_identity2) + "\n")  
    fob2.write("low percent identity cutoff: " + str(low_identity2) + "\n")
    fob2.write("high length identity cutoff: " + str(high_length2) + "\n")
    fob2.write("low length identity cutoff: " + str(low_length2) + "\n")
    fob2.write("write query file name and subject IDs: " + w_names + "\n")
    fob2.write("write query and subject lengths, and percent coverage: " + w_lengths + "\n")
    fob2.write("write number of identical residues and percent identity: " + w_identities + "\n")
    fob2.write("write alignment: " + w_alignment + "\n")
    fob2.write("alignment output length: " + str(output_alignment_length) + "\n")
    fob2.close()
    
def write_header():
    fob.write("~~~~~~~~~~~~~~~~~~~~~~~~~~ Output ~~~~~~~~~~~~~~~~~~~~~~~~~" + "\n")
    fob.write("Result " + filename + "-" + str(counter).zfill(3) + "\n")
    fob.write("\n")

def write_names():
    fob.write("query file: " + filename + "\n")
    fob.write("\n")
    fob.write("sequence: " + "\n")
#    initial_alignment_length()           ###Weird
    first = 0                             ###Weird
    second = output_alignment_length      ###Weird
    fob.write(alignment.title[first:second] + "\n") #Comment out to write full title
    
#    for i in range((int(len(alignment.title) // output_alignment_length) + 1)): #Uncomment for full title
#        fob.write(alignment.title[first:second] + "\n")    #Uncomment for full title
#        first = first + output_alignment_length            #Uncomment for full title
#        second = second + output_alignment_length          #Uncomment for full title
    fob.write("\n") 

def write_lengths():
    coverage = format(float(alignment.length) / float(query_length) * 100.0, '.2f')
    fob.write("query length: " + str(query_length) + "\n")     
    fob.write("alighment length: " + str(alignment.length) + "\n")
    fob.write("percent query coverage: " + str(coverage) + "%" + "\n")
    fob.write("\n")

def write_identities():
    percent_identity2 = format(percent_identity * 100, '.2f')
    fob.write("identical residues: " + str(identical_residues) + "\n")
    fob.write("percent identity: " + str(percent_identity2) + "%"+ "\n")
    fob.write("\n")
    
def write_alignment_markers():
    """This makes the '1234567890' string above the alignment
    it uses a couple of loops so it may change size according
    to what the user enters as an alignment length"""
    marker_out = ""    
    for i in range(output_alignment_length // 10):
        marker_out = marker_out + "1234567890"
        
    count = 0
    for i in range(output_alignment_length % 10):
        count = count +1
        marker_out = marker_out + str(count)
        
    fob.write(marker_out + "\n")
        
def write_alignment():
    write_alignment_markers()
#    initial_alignment_length()       ###Weird
    first = 0                         ###Weird
    second = output_alignment_length  ###Weird
    for i in range((int(len(hsp.match) // output_alignment_length) + 1)):
        fob.write(hsp.query[first:second] + "\n")
        fob.write(hsp.match[first:second] + "\n")
        fob.write(hsp.sbjct[first:second] + "\n")
        fob.write("\n")
        
        first = first + output_alignment_length
        second = second + output_alignment_length

def write_stats():
    fob3.write(filename + " gave " + str(counter) + " results" + "\n")

def write_stats2():
    fob3.write("Total alignments: " + str(total) + "\n")


#### program start #### 
                                   #### read parameters from parse.inp
read_param_file()

#check_output_dir()
                                   #### write parameters to a .def file based on output root
#write_param_file()
                                   #### open output file to write out the results and stats
#init_alignment_length()

#fob = open(output_file_path + "/" + output_file_root + ".out", 'w')  # Comment/Uncomment - write a single output file, named after the output file root
fob3 = open(output_file_root + ".stats", 'w')
                                   #### gets the length of query and stores to a variable
total = 0
for filename in glob.glob(input_file_names):
    record = SeqIO.read(filename, "fasta")
    query_length = len(record)
                                   #### compare hits to current input query
                                   #### define the handle
    filename2 = filename + ".xml" 
    result_handle = open(filename2)
    
    blast_record = NCBIXML.read(result_handle)
                                   #### write query file name
    counter = 0

                                   #### screen blast output records against parameters
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            alignment_length = alignment.length
            identical_residues = hsp.identities
            percent_identity = float(identical_residues) / float(query_length) 

            cond1 = percent_identity <= high_identity 
            cond2 = percent_identity > low_identity
            cond3 = alignment_length <= query_length * high_length
            cond4 = alignment_length > query_length * low_length
                                   #### write blast output that passes screens
            if cond1 and cond2 and cond3 and cond4:
                counter += 1
#                fob = open(output_file_path + "/" + filename + "-" + str(counter).zfill(3) + ".out", 'w') # Comment/Uncomment - write individual output files, with filenames based on input files.
                
#                 if w_header[:1] == "Y" or w_names[:1] == "y":
#                     write_header()
#                 if w_names[:1] == "Y" or w_names[:1] == "y":
#                     write_names()
#                 if w_lengths[:1] == "Y" or w_lengths[:1] == "y":
#                     write_lengths()
#                 if w_identities[:1] == "Y" or w_identities[:1] == "y":
#                     write_identities()
#                 if w_alignment[:1] == "Y" or w_alignment[:1] == "y":
#                     write_alignment()                
    write_stats()
    total = total + counter
write_stats2()
#fob.close()
fob3.close()

end=time.time()
print end - start