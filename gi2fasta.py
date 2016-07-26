import sys
#from Bio import Entrez
#from Bio import SeqIO

user_email = "" # User must supply email here to access NCBI api
# Add error message in the event no email address is supplied 
if user_email == "":
    sys.exit("Error: Please supply your email address to line 5 of gi2fasta.py")

Entrez.email = user_email

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

