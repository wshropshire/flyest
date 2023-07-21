#!/bin/bash

# Check if the required arguments are provided
if [ $# -ne 3 ]; then
  echo "Usage: $0 <flye_directory> <min_coverage_depth> <min_length>"
  exit 1
fi

# Assign the command-line arguments to variables
flye_directory=$1
min_coverage_depth=$2
min_contig_length=$3

# Check if the flye directory exists
if [ ! -d "$flye_directory" ]; then
  echo "Error: Flye directory does not exist."
  exit 1
fi

# Read the assembly_info.txt file
assembly_info_file="$flye_directory/assembly_info.txt"
if [ ! -f "$assembly_info_file" ]; then
  echo "Error: assembly_info.txt file not found in the flye directory."
  exit 1
fi

# Extract sequence names with coverage depth greater than or equal to the minimum coverage depth
sequence_names=$(awk -v min_cov="$min_coverage_depth" -v min_length="$min_contig_length" '(NR > 1 && $3 >= min_cov && $2 >= min_length) { print $1 }' "$assembly_info_file")
discard_names=$(awk -v min_cov="$min_coverage_depth" -v min_length="$min_contig_length" '(NR > 1 && $3 < min_cov || $2 < min_length) { print $1 }' "$assembly_info_file")

# Check if any sequences meet the coverage depth requirement
if [ -z "$sequence_names" ]; then
  echo "No sequences found with coverage depth greater than or equal to $min_coverage_depth."
  exit 0
fi

# Read the draft assembly file
draft_assembly_file="$flye_directory/assembly.fasta"
if [ ! -f "$draft_assembly_file" ]; then
  echo "Error: assembly.fasta file not found in the flye directory."
  exit 1
fi

# Create the filtered fasta file
filtered_fasta_file="$flye_directory/filtered.fasta"
discarded_fasta_file="$flye_directory/discarded.fasta"

# Process the filtered sequences
if [ -n "$sequence_names" ]; then
  awk -v seq_names="$sequence_names" '
    BEGIN {
      split(seq_names, names, " ")
      flag = 0
      seq = ""
    }
    /^>/ {
      if (flag && seq != "") {
        print seq
        seq = ""
      }
      flag = ($1 ~ names[1])
      for (i = 2; i <= length(names); i++) {
        if ($1 ~ names[i]) {
          flag = 1
          break
        }
      }
      if (flag) {
        print
      }
    }
    flag {
      if ($0 !~ /^>/) {
        seq = seq $0
      }
    }
    END {
      if (flag && seq != "") {
        print seq
      }
    }
  ' "$draft_assembly_file" > "$filtered_fasta_file"
  echo "Filtered fasta file created: $filtered_fasta_file"
fi

# Process the discarded sequences
if [ -n "$discard_names" ]; then
  awk -v discard_names="$discard_names" '
    BEGIN {
      split(discard_names, names, " ")
      flag = 0
      seq = ""
    }
    /^>/ {
      if (flag && seq != "") {
        print seq
        seq = ""
      }
      flag = ($1 ~ names[1])
      for (i = 2; i <= length(names); i++) {
        if ($1 ~ names[i]) {
          flag = 1
          break
        }
      }
      if (flag) {
        print
      }
    }
    flag {
      if ($0 !~ /^>/) {
        seq = seq $0
      }
    }
    END {
      if (flag && seq != "") {
        print seq
      }
    }
  ' "$draft_assembly_file" > "$discarded_fasta_file"
  echo "Discarded fasta file created: $discarded_fasta_file"
else
  echo "No sequences found below the minimum coverage depth or minimum contig length."
fi

