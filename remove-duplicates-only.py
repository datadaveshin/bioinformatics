"""
Dave Shin
2015-07-22 5:36:51 PM
Description: Removes duplicate fasta files, exports list of removed duplicates and number of cysteines
Modified from: remove-duplicates3D.py by David Shin and Derek Green
Usage: python remove-duplicates.py <input file> <output file>
Example: python remove-duplicates.py xac-fasta.out xac-fasta
Version changes:
Removed counting cysteines and glycosylation sites
Now just removes duplicate entries, however, it seems Blast avoids these, this is just for check
"""
import sys
                                        ## read in input file and output file root names
infilename = sys.argv[1]  
outroot = sys.argv[2]

# 
proteins = {}
file = open(infilename)
outfile = open(outroot + '.kept', 'w')
removedfile = open(outroot + '.removed','w')

for line in file:
    if line.startswith('>'):
        header = line.strip()
    else:
        seq = line.strip()
        print seq
        if seq in proteins:
            print "excluded the following duplicate file:", header
            removedfile.write("excluded the following duplicate file: " + header + '\n')
        else:
            proteins[seq] = header
            outfile.write(header + "\n")
            outfile.write(seq + "\n")
outfile.close()
removedfile.close()
