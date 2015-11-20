# Bioinformatics
Repository of shell and Python scripts for bioinformatics.

## deltablast_xml2aln.csh
Takes .xml deltablast output from NCBI's site, fetches full-length hits in fasta format according to user defined % identity and protein length requirements. The fasta files are screened to remove any accidental duplicates, then are aligned using clustalo.

Dependencies/Requirements: 
1) Python 2.7
2) The Biopython module
3) A local copy of clustal O on your computer
