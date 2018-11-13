# Bioinformatics
###### Repository of Shell and Python scripts for bioinformatics.

## Retrieve Full-Length .fasta Files following Deltablast and Generate a Clustal Alignment
### deltablast_xml2aln.csh 
#### Description:
`deltablast_xml2aln.csh` takes `.xml` deltablast output from NCBI's site, fetches full-length hits in fasta format according to user defined % identity and protein length requirements. The fasta files are screened to remove any accidental duplicates, then are aligned using `clustalo`.
#### Requirements: 
1. `Python 2.7`
2. The `Biopython` Python module - http://biopython.org/wiki/Download or easier through an anaconda/conda setup from https://www.continuum.io/downloads
3. A local copy of clustal omega - `clustalo` - on your computer - available from http://www.clustal.org/omega/

#### References:
1. Cock PA, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Friedberg I, Hamelryck T, Kauff F, Wilczynski B, de Hoon MJL (2009) Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25, 1422-1423
1. Sievers F, Wilm A, Dineen DG, Gibson TJ, Karplus K, Li W, Lopez R, McWilliam H, Remmert M, Söding J, Thompson JD, Higgins DG (2011) Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega. Molecular Systems Biology 7:539 doi:10.1038/msb.2011.75

## Protein Motif, Site Finder and Statistics from an Alignment
### pattern_finder.py
#### Description:
`pattern_finder.py` takes your clustal alignment and a motif pattern to search for as input and produces a new alignment with several style options showing where your motif pattern exists and does not exist relative to the parent sequence. Statistics are given after printing the alignment. 

#### Requirements: 
1. `Python 2.7`
1. A clustal omega alignment
