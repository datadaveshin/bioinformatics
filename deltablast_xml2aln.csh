#!/bin/csh
# Copyright (c) 2015-11-19 7:20:47 PM, David Shin 
# Description: 
# Modified from: 
# Usage: 
# Example: ./launch parse-95-60-115-85.inp  pepck
#################################################

# Trim the blast xml output using the supplied parse.inp file
# Code from parse-launch2.com

python deltablast_parse.py $1 $2

# Print out how many alignments after parsing
# Code from parse-stats.com

python deltablast_parse_stats.py $1

# Change directories into the results directory

cd $2-*_results

# Gives the gi numbers in a single list so can use entrez to get the fasta files
# Code from launch-parse-select-keep-all.com

python parse_select_keep_all.py $2

# Gets fasta files from gi numbers in the ###-gi.out list the results directory
# Note, this also gets a description of the file
# Code from launch-gi-to-fasta.com 

python gi2fasta.py $2-gi.out $2-fasta.out

# Checks for and removes duplicate entries based on sequence
# May not be needed as it appears blast doesn't return duplicate sequences
# Code from launch-remove-duplicates.com but updated the remove-duplicates python script
# To remove searching for cysteines and glycosylation sites, that should be separate script

python remove_duplicates_only.py $2-fasta.out $2-fasta

# Add the original search fasta file to the top of the "kept" list

cat ../$2 $2-fasta.kept > $2-fasta2.out

# Make clustal alignments
# Code from launch-clustalo.com

clustalo -i $2-fasta2.out --outfmt clu --residuenumber --output-order=input-order --wrap 10000 -v --outfile $2-fasta.aln

