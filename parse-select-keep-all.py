"""
Dave Shin
Thu Dec 12 02:20:54 PST 2013
Description: Used to select alignments coming out of deltablast-parse.py \n
Modified from: 
Usage: python parse-select.py <deltablast-parse output root>
    Then hit return to save an alignment and associated gi code,
    or any other key followed by return to reject an alignment.
Example: python parse-select.py xaa
Version 2m
2m: updated header
2l: now <return> only will let you keep an alignment, any other key followed by return is a rejection
"""

import glob                     # used to read in files from a directory
import os                       # used to clear screen
import sys                      # used to read in file names from command line
                                        ## slices title to get gi number
def find_between( s, first, last ):
    try:
        start = s.index( first )+ len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return""

                                        ## writes out just the gi number after going thru find_between()
def write_gi():
    from_file = open(filename, 'r')
                                    # reads the file into a list
    my_file_list = []
    my_file_list = from_file.readlines()
                                   
    gi = my_file_list[6]            # gets the gi number from line 6
    gi2 = find_between(gi, "|", "|")
    fob.write(gi2 + "\n")
    from_file.close()

                                    #### program start

infile_root = sys.argv[-1]

fob = open(str(infile_root) + '-gi.out', 'w')

input_file_names = str(infile_root + "-???.out")                            

for filename in glob.glob(input_file_names):
    os.system("clear")
    f = open(filename, 'r')
    print f.read()
    
#    decision = raw_input("Keep this sequence? [Return = yes, any key and return = no]")
#    if decision == "":
    write_gi()
    f.close()

fob.close()
