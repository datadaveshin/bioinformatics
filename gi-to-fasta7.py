import sys
from Bio import Entrez
from Bio import SeqIO
Entrez.email = "davidsshin@lbl.gov"

infilename = sys.argv[1]
outfilename = sys.argv[2]

with open(infilename) as f:
    gi_numbers=', '.join(line.rstrip() for line in f)

handle = Entrez.efetch(db="protein", rettype="fasta", retmode="text", id=gi_numbers)
records = SeqIO.parse(handle, "fasta")

fout = open(outfilename, 'w') 

for record in records:
    #print ">" + record.seq
    #print record.id
    print record.description
    #print record.seq	
    fout.write(">" + str(record.description) + "\n")
    fout.write(str(record.seq) + "\n")
fout.close()
#for seq_record in SeqIO.parse(record, "fasta"):
#    print seq_record.id

#fob2.write("high percent identity cutoff: " + str(high_identity2) + "\n")

