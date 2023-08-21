#!/usr/bin/env python3

from Bio import SeqIO
import sys
import os

# Check if the required arguments are provided
if len(sys.argv) != 4:
    print("Usage: {} <flye_directory> <min_coverage_depth> <min_length>".format(sys.argv[0]))
    sys.exit(1)

# Assign command-line arguments to variables
flye_directory = sys.argv[1]
min_coverage_depth = float(sys.argv[2])
min_contig_length = int(sys.argv[3])

# ... (other initial checks)

# Read the assembly_info.txt file and extract sequence names
sequence_names = []
with open(os.path.join(flye_directory, "assembly_info.txt"), "r") as info_file:
    for line in info_file:
        if not line.startswith("#"):
            fields = line.strip().split("\t")
            seq_name = fields[0]
            cov = float(fields[2])
            length = int(fields[1])
            if cov >= min_coverage_depth and length >= min_contig_length:
                sequence_names.append(seq_name)

# Read the draft assembly file
draft_assembly_file = os.path.join(flye_directory, "assembly.fasta")
if not os.path.isfile(draft_assembly_file):
    print("Error: assembly.fasta file not found in the flye directory.")
    sys.exit(1)

# Create the filtered and discarded fasta files
filtered_fasta_file = os.path.join(flye_directory, "filtered.fasta")
discarded_fasta_file = os.path.join(flye_directory, "discarded.fasta")

# Process the sequences
filtered_seqs = []
discarded_seqs = []

with open(draft_assembly_file, "r") as fasta_file:
    current_seq = []
    for line in fasta_file:
        if line.startswith(">"):
            if current_seq:
                if current_seq[0][1:] in sequence_names:
                    filtered_seqs.append(current_seq)
                else:
                    discarded_seqs.append(current_seq)
            current_seq = [line.strip()]
        else:
            current_seq.append(line.strip())

# Process the last sequence
if current_seq:
    if current_seq[0][1:] in sequence_names:
        filtered_seqs.append(current_seq)
    else:
        discarded_seqs.append(current_seq)

# Write filtered sequences to the output file
with open(filtered_fasta_file, "w") as f:
    for seq in filtered_seqs:
        f.write("\n".join(seq) + "\n")

print("Filtered fasta file created:", filtered_fasta_file)

# Write discarded sequences to the output file
with open(discarded_fasta_file, "w") as f:
    for seq in discarded_seqs:
        f.write("\n".join(seq) + "\n")

print("Discarded fasta file created:", discarded_fasta_file)

