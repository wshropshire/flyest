#!/bin/bash

# Usage output
usage() {

    echo "Usage: $0 [-i <input_file>] [-o <output_dir>] [-p <prefix>]"
    echo "       -i <input_file>: Polished consensus genome assembly"
    echo "       -o <output_dir>: Output directory (default: 'genome_analysis')"
    echo "       -p <prefix>: Prefix for the output files (default: 'sample')"
	echo "STILL EXPERIMENTAL: DOESN'T WORK DUE TO DEPENDENCY CONFLICTS"
}

# Default values
output_dir="genome_analysis"
prefix="sample"

# Parse input options
while getopts ":i:o:p:" opt; do
    case $opt in
        i)
            input_file=$OPTARG
            ;;
        o)
            output_dir=$OPTARG
            ;;
        p)
            prefix=$OPTARG
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            usage
            exit 1
            ;;
    esac
done

# Check if the input file is provided
if [ -z "$input_file" ]; then
    usage
    exit 1
fi

# Check if the output directory already exists
if [ -d "$output_dir" ]; then
    echo "Error: Output directory '$output_dir' already exists. Please choose a different output directory."
    exit 0
else
	echo ""
        echo "Making directories and log files"
	mkdir $output_dir
	cd "$output_dir" || exit
        # Logging
        LOG_DIR=$OUT_DIR/logs
        mkdir $LOG_DIR
        LOG_NAME="$LOG_DIR/${SAMPLE}_flyer_$(date +"%Y-%m-%d").txt"
        echo -e "${SAMPLE}_flyer_log_$(date +"%Y-%m-%d")\n" >> $LOG_NAME
        exec &> >(tee -a "$LOG_NAME")
	#### Program versions
        echo -e "### Program versions:\n" >> $LOG_NAME
        echo "Flye:"$(flye --version) >> $LOG_NAME


# MLST
mlst_dir="mlst"
mkdir -p "$mlst_dir"
cd "$mlst_dir" || exit
echo "Running MLST..."
mlst "$input_file" > "$prefix"_mlst.txt

# BUSCO QC
busco_dir="busco"
mkdir -p "$busco_dir"
cd "$busco_dir" || exit
echo "Running BUSCO for quality control..."
busco -i "../$input_file" -o "$prefix"_busco_output -l bacteria_odb10 -m genome

# Prokka annotation
cd ..
prokka_dir="prokka_annotation"
mkdir -p "$prokka_dir"
cd "$prokka_dir" || exit
echo "Running Prokka annotation..."
prokka --compliant --outdir prokka_output --prefix "$prefix"_prokka "$input_file"

# AMRFinderPlus
cd ..
amrfinderplus_dir="amrfinderplus"
mkdir -p "$amrfinderplus_dir"
cd "$amrfinderplus_dir" || exit
echo "Running AMRFinderPlus..."
amrfinderplus -n "../$input_file" -o "$prefix"_amrfinderplus_output

echo "Analysis completed."

