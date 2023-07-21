#!/bin/bash

#Set the environment variables
export FASTA_FILE="$1"
export FASTQ_FILE="$2"
export OUTPUT_DIR="$3"
export SAMPLE="$4"
export THREADS="$5"
export MINIMAP2=minimap2
export SAMTOOLS=samtools
export BEDTOOLS=bedtools

# Extract contig names and sequences from the multi-fasta file
awk '/^>/{if (seq) { print seq; } printf("%s\t",$0); seq=""; next; } { seq = seq $0 } END { print seq; }' "$FASTA_FILE" > $OUTPUT_DIR/contigs.txt

# Calculate the base pair count for each contig
while IFS=$'\t' read -r header sequence; do
    contig=$(echo "$header" | sed 's/^>//')
    bp_count=$(echo "$sequence" | tr -d '\n' | wc -c)
    echo "Contig: $contig, Base Pairs: $((bp_count))"
done < $OUTPUT_DIR/contigs.txt > $OUTPUT_DIR/${SAMPLE}_qc.txt

# Complete alignment using minimap2 
bam_file=$OUTPUT_DIR/${SAMPLE}_consensus_ONT.bam
$SAMTOOLS faidx $FASTA_FILE
$MINIMAP2 -t $THREADS -ax map-ont $FASTA_FILE $FASTQ_FILE | $SAMTOOLS 'sort' -@ $THREADS > $bam_file
$SAMTOOLS 'index' $bam_file

$BEDTOOLS genomecov -ibam $bam_file -d > $OUTPUT_DIR/${SAMPLE}_coverage.bed


# Extract the coverage depths to a separate file
awk '{ print $3 }' $OUTPUT_DIR/${SAMPLE}_coverage.bed > $OUTPUT_DIR/coverage_depths.txt

# 4: Calculate the mean coverage depth using awk
mean=$(awk '{ sum += $1 } END { print sum / NR }' $OUTPUT_DIR/coverage_depths.txt)
echo "" >> $OUTPUT_DIR/${SAMPLE}_qc.txt
echo "Mean Coverage Depth: $mean" >> $OUTPUT_DIR/${SAMPLE}_qc.txt

# Step 5: Calculate the standard deviation using awk
stddev=$(awk -v mean="$mean" '{ sum += ($1 - mean) * ($1 - mean) } END { print sqrt(sum / NR) }' $OUTPUT_DIR/coverage_depths.txt)
echo "Standard Deviation: $stddev" >> $OUTPUT_DIR/${SAMPLE}_qc.txt
# Clean up temporary files
rm $OUTPUT_DIR/contigs.txt $OUTPUT_DIR/coverage_depths.txt
