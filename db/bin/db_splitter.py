#!/usr/bin/env python

# This script allows to split blast databases such as nr or nt
# based into separate databases based on the taxonomic categories
# Requires the output of gen_db_headers for the blastdb

# Assumes the input is the fasta output fromm the blastdbcmd with the
# -ctrl_a option

import sys
import re
from collections import deque


def main(): 
    db_in = sys.stdin
    db_out = sys.stdout
    
    db_name = "nr"
    
    seq_label_file = "../../db_headers/nr_headers"
    seq_labels_fh = open(seq_label_file)
    
    
    #NCBI taxonomy db division categories are 12
    # BCT, INV, MAM, PHG, PLN, PRI, ROD,SYN, UNA, VRL, VRT, ENV
    #output_categories = ["MAM", "BCT", "VRL", "OTR"]
    
    db_to_output_category_map = {   "MAM": "MAM", "ROD":"MAM", "PRI":"MAM",
                                    "BCT": "BCT",
                                    "PHG": "VRL", "VRL":"VRL",
                                    "INV": "OTR","PLN":"OTR", "SYN": "OTR", "UNA":"OTR", "VRT":"OTR", "ENV": "OTR" 
                                 }
    #PARAMS
    label_to_keep = "MAM"
    
    labels_in_memory = 10000 #how many labels to bulk load from the header file
    fasta_lines_in_memory = 10000
    
    
    #Load 
    nSeqs_processed  = 0
    include_sequence = False
    
    label_reader = bulkLoadSequenceLabels(seq_labels_fh,labels_in_memory)
    seq_labels = label_reader.next()
    
    for fasta_lines in bulkLoadFasta(db_in, fasta_lines_in_memory):
        for fasta_line in fasta_lines:
            #If fasta sequence header line
            if re.match(fasta_line,"^>"):
                include_sequence = belongs_to_category(fasta_line, label_to_keep, seq_labels.popleft() )
                nSeqs_processed += 1
                if include_sequence:
                    db_out.write(fasta_line)
                
                if nSeqs_processed % labels_in_memory == 0: #If loaded labels are empty, load next batch of labels
                    seq_labels = label_reader.next()
            else: #If it is nucleotide/aminoacid sequnce
                if include_sequence:
                    db_out.write(fasta_line)
        
def belong_to_category(fasta_header, label_to_keep,  seq_label ):
    #Verify if label corresponds to the read header:
    if fasta_header.find(seq_label[0]) > -1:
        return label_to_keep == seq_label[1]
    else:
        sys.stderr("ERROR! Read label does not correspond to fasta line")
        raise Exception("Fail :(")

#Lazy loading methods!
def bulkLoadSequenceLabels(label_file ,nLines = 10000):
    while True:
        data = deque()
        for i in range(nLines):
            line = dbfile.readline()
            if not line:
                break
            #Extract gi and tax. category (cols 0 and 2)
            fields = line.rstrip("\n").split("\t")
            data.append( (fields[0],fields[2] ) )
        if not data:
            break
        yield data
        

def bulkLoadFasta(db_file,nLines= 10000):
    while True:
        data = []
        for i in range(nLines):
            line = dbfile.readline()
            if not line:
                break
            data.append(line)
        if not data:
            break
        yield data
        
if __name__ == '__main__':
    main()