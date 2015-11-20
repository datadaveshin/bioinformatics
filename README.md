# Bioinformatics
###### Repository of shell and Python scripts for bioinformatics.

## Deltablast hits to full-length .fasta files to clustal alignment
deltablast_xml2aln.csh takes .xml deltablast output from NCBI's site, fetches full-length hits in fasta format according to user defined % identity and protein length requirements. The fasta files are screened to remove any accidental duplicates, then are aligned using clustalo.

##### Dependencies/Requirements: 
1. Python 2.7
2. The Biopython module
3. A local copy of clustal O on your computer

## Pattern Finder 
pattern_finder.py takes your clustal alignment and a motif pattern to search for and produces a new alignment with several style options showing where your motif pattern exists and does not exist relative to the parent sequence. Statistics are given after printing the alignment. 

##### Dependencies/Requirements: 
1. Python 2.7
2. The Biopython module
3. A local copy of clustal O on your computer
