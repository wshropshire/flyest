#!/bin/bash

### Description ----------------------------------------------------------------

USAGE="
-- Nanopore + Illumina hybrid consensus polishing pipeline -- flyest

usage: $(basename "$0" .sh) [-h] [--nano-corr] [-n value] [-l value] [--meta] [-m string] [--nid value] [--nlen value] [--no-qc] [-t value] [-s string]
(-i string -1 string -2 string -o string)

where:
	-i|--input	File containing long-read fastq reads [Can be gzipped] --REQUIRED.
	-1|--pe1	Paired-end forward short-reads --REQUIRED.
	-2|--pe2	Paired-end reverse short-reads --REQUIRED.
	--nano-corr	Long-read input quality [Default=--nano-hq for ONT high-quality reads; Only use --nano-corr for pre Guppy-v5.x basecalled reads].
	-n|--min	Minimum coverage depth per sequence as a fraction of the overall draft assembly mean coverage depth [Default=0.2; i.e. 20%]. Use a value of '0.001' for no minimum coverage depth.
	-l|--len	Minimum contig length per sequence [Default=1000].
	--meta	Meta option for Flye [Default is not set] - Good for uneven/low coverage assemblies.
	-m|--mod	Medaka model [Default=r1041_e82_400bps_sup_v4.2.0].
	--nid	Nucmer min_id parameter for removing highly similar contigs based on percent minimum nucleotide identity shared [Default = 95].
	--nlen	Nucmer min_length parameter for removing highly similar contigs based on percent minimum contig length shared [Default = 90].
	--no-qc	Disable quick QC script [Default is on].
	-o|--outdir	Output directory --REQUIRED.
	-t|--threads	Number of threads to use [Default=1].
	-s|--sample	Sample prefix [DEFAULT=SAMPLE].
	--keep-tmp      Preserve all intermediate files within pathway for troubleshooting [Default is remove] - CURRENTLY UNENABLED.
	-h|--help|-u|--usage	Show this help text.
"

### Terminal Arguments ---------------------------------------------------------

if [[ $# -eq 0 ]];then
	echo "$USAGE"
	exit 1
fi

# Conditional default environment variables
EXECUTE_QUICK_QC=true
FLYE_READ_MOD="--nano-hq"
META=""
MIN_COV='0.2'
MIN_LENGTH='1000'
MEDAKA_MODEL="r1041_e82_400bps_sup_v4.2.0"
NUC_ID="95"
NUC_LEN="90"
THREADS=1
SAMPLE="SAMPLE"
KEEP_TMP=false
VERSION="1.0.0"

while [[ $# > 0 ]]
do
	key="$1"

	case $key in
		-i|--input)
			export FASTQ_FILE="$2"
			shift
			;;
		-1|--pe1)
			export PE1="$2"
			shift
			;;
		-2|--pe2)
			export PE2="$2"
			shift
			;;
		--nano-corr)
			FLYE_READ_MOD="--nano-corr"  
			;;
		-n|--min)
			export MIN_COV="$2"
			shift
			;;
		-l|--len)
			export MIN_LENGTH="$2"
			shift
			;;
		--meta)
			META="--meta"  
			;;
		-m|--mod)
			export MEDAKA_MODEL="$2"
			shift
			;;
		--nid)
			export NUC_ID="$2"
			shift
			;;
		--nlen)
			export NUC_LEN="$2"
			shift
			;;
		--no-qc)
            		EXECUTE_QUICK_QC=false
            		;;
		-o|--outdir)
			export OUT_DIR="$2"
			shift
			;;
		-t|--threads)
			export THREADS="$2"
			shift
			;;
		-s|--sample)
			export SAMPLE="$2"
			shift
			;;
		--keep-tmp)
			KEEP_TMP=true
			shift
			;;
		-h|--help|-u|--usage)
			echo "$USAGE"  
			exit 0
			;;
		*)
			echo "Unknown option $1"
			exit 1
			;;
	esac
	shift
done

# Handle Missingness for requried variables
MISSING="is missing but required. Exiting."
if [ -z ${FASTQ_FILE+x} ]; then echo "-i $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${PE1+x} ]; then echo "-1 ${MISSING}. If creating long-read only assembly, use flyer script"; echo "$USAGE"; exit 1; fi;
if [ -z ${PE2+x} ]; then echo "-2 ${MISSING}. If creating long-read only assembly, use flyer script"; echo "$USAGE"; exit 1; fi;
if [ -z ${OUT_DIR+x} ]; then echo "-o $MISSING"; echo "$USAGE"; exit 1; fi;

### Create directories and logging ---------------------------------------------------------

# Record the start time
START_TIME=$(date +%s)

# Set script directory pathway variable
SCRIPT_DIR="$(dirname "$0")"

# Create directory if not currently existing and proceed to logging. If directory existing, echo usage
if [ -d "$OUT_DIR" ]
        then echo "$OUT_DIR exists. Exiting..."
        exit 0
else
	echo ""
        echo "Making directories and log files"
	mkdir $OUT_DIR
        # Logging
        LOG_DIR=$OUT_DIR/logs
        mkdir $LOG_DIR
        LOG_NAME="$LOG_DIR/${SAMPLE}_flyest_$(date +"%Y-%m-%d").txt"
        echo -e "${SAMPLE}_flyest_log_$(date +"%Y-%m-%d")\n" >> $LOG_NAME
        exec &> >(tee -a "${LOG_NAME}")
        #### Program versions
        echo "### Program versions:" >> $LOG_NAME
        echo "Flyest:"$VERSION >> $LOG_NAME
	echo "Flye:"$(flye --version) >> $LOG_NAME
        echo "Medaka:"$(medaka --version | cut -d' ' -f2) >> $LOG_NAME
        echo "dnaapler:"$(dnaapler --version | cut -d' ' -f3) >> $LOG_NAME
	echo "berokka:"$(berokka --version | cut -d' ' -f2) >> $LOG_NAME
        echo "pilon:"$(pilon --version | cut -d' ' -f3) >> $LOG_NAME 
	echo "polypolish:"$(polypolish --version | cut -d' ' -f2) >> $LOG_NAME
        echo -e "polca:Derived from MaSuRCA-v4.1.0\n" >> $LOG_NAME
	echo "### Settings:"
        echo "Input long-reads: $FASTQ_FILE"
	echo "Input short-reads: $PE1 and $PE2"
	echo "Flye read quality: $FLYE_READ_MOD"
	if [[ "${META}" != "" ]]; then
                echo "Flye meta option enabled: Yes"
        else
                echo "Flye meta option enabled: No"
        fi
        echo "Minimum contig coverage depth threshold as fraction of assembly mean coverage depth: $MIN_COV"
        echo "Minimum contig length: $MIN_LENGTH"
	echo "Medaka model: ${MEDAKA_MODEL}"
        echo "Output directory: $OUT_DIR"
        echo "Prefix: $SAMPLE"
        echo "Threads: ${THREADS}"
	if [[ "$EXECUTE_QUICK_QC" = false ]]; then
		echo -e "Quick QC enabled: No\n"
        else
                echo -e "Quick QC enabled: Yes\n"
        fi
fi

### Pipeline -----------------------------------------------------------

#### Set executable file environment variables - Need to either create virtual environment or replace commands with absolute pathways/configurements - highly recommend the former given that tools such as Medaka require further dependencies not listed here
export FLYE=flye
export MEDAKA=medaka_consensus
export PILON=pilon
export MINIMAP2=minimap2
export SAMTOOLS=samtools
export BEROKKA=berokka
export BWA=bwa
export POLY_INSERT=polypolish_insert_filter.py
export POLYPOLISH=polypolish
export POLCA=polca_mod.sh 

#### Draft Assembly
echo -e "Step 1: Creating draft assembly with Flye\n"

$FLYE \
	$FLYE_READ_MOD $FASTQ_FILE \
	--threads $THREADS \
	--out-dir $OUT_DIR/${SAMPLE}_flye_assembly \
	$META
echo ""
echo -e "Flye draft assembly complete\n"

# Calculate the minimum coverage depth based on default or user specified fraction of the mean coverage draph assembly
MEAN_COV=$(grep "Mean coverage:" $OUT_DIR/${SAMPLE}_flye_assembly/flye.log | cut -d$'\t' -f3)
echo -e "Calculating the minimum coverage depth based on ${MIN_COV} fraction of draft assembly coverage depth=${MEAN_COV}\n"

# Use command line calculator 'bc' to do arithmetic with float variables
MIN_COV_FLOAT=$(echo "$MIN_COV" | bc -l)
CLEAN_COV_INPUT=$(echo "$MEAN_COV * $MIN_COV_FLOAT" | bc -l)
echo -e "Minimum coverage depth will be ${CLEAN_COV_INPUT}\n"

# Remove sequences in draft assembly with coverage depth less than minimum: Default = 20% of mean coverage depth of draft assembly
echo -e "Step 2: Evaluating minimum coverage depth and minimum contig length per sequence in draft assembly\n"
python $SCRIPT_DIR/flye_draft_clean.py $OUT_DIR/${SAMPLE}_flye_assembly $CLEAN_COV_INPUT $MIN_LENGTH

# TODO: Consider dump polishing if contig numbers exceed a number

#Set the paths and names for input and output files for dnaapler1
FASTA_INPUT1=$OUT_DIR/${SAMPLE}_flye_assembly/filtered.fasta
FASTA_OUTPUT1=$OUT_DIR/${SAMPLE}_"dnaapler1.fasta"
FUNCTION1="mystery"

echo -e "First rotation of assembly using random gene with dnaapler\n"
# Call dnaapler mystery for first rotation of draft assembly 
source $SCRIPT_DIR/dnaapler_script.sh "$FASTA_INPUT1" "$FASTA_OUTPUT1" "$FUNCTION1" "$OUT_DIR"

echo ""
echo -e "Step 3: Running Medaka polish using $MEDAKA_MODEL as a model\n"
$MEDAKA \
	-i $FASTQ_FILE -d $FASTA_OUTPUT1 -o $OUT_DIR/${SAMPLE}_medaka -m $MEDAKA_MODEL -t $THREADS

MEDAKA_OUTPUT=$OUT_DIR/${SAMPLE}_medaka/consensus.fasta

echo ""
echo -e "Step 4: Run Berokka and modified clean script to trim overlaps from draft assembly and remove identical contigs based on the following parameters: numcer min_id = $NUC_ID and nucmer min_length = ${NUC_LEN}\n"

# Run berokka
$BEROKKA \
        --outdir $OUT_DIR/${SAMPLE}_berokka_clean $MEDAKA_OUTPUT

BEROKKA_OUTPUT=$OUT_DIR/${SAMPLE}_berokka_clean/02.trimmed.fa

# Run modified script adapted from circlator clean (Hunt et al. 2015 Genome Biology)
python $SCRIPT_DIR/clean.py --input $BEROKKA_OUTPUT --outdir $OUT_DIR/${SAMPLE}_berokka_clean --prefix ${SAMPLE}_clean --min_length $MIN_LENGTH --nucmer_min_id $NUC_ID --nucmer_min_length $NUC_LEN

# Set the paths and names for input and output files for dnaapler2 - include bam file for pilon
FASTA_INPUT2=$OUT_DIR/${SAMPLE}_berokka_clean/${SAMPLE}_clean.fasta
FASTA_OUTPUT2=$OUT_DIR/${SAMPLE}_"dnaapler2.fasta"
BAM_FILE=$OUT_DIR/${SAMPLE}_dnaapler2.bam

# Call dnaapler mystery for second rotation of draft assembly
echo ""
echo -e "Second rotation of assembly using random gene with dnaapler\n"
$SCRIPT_DIR/dnaapler_script.sh "$FASTA_INPUT2" "$FASTA_OUTPUT2" "$FUNCTION1" "$OUT_DIR"

# Call pilon for final long-read polish - Note that it is a memory hog 
echo ""
echo -e "Step 5: Performing final long-read polish with pilon\n"
$MINIMAP2 -t $THREADS -ax map-ont $FASTA_OUTPUT2 $FASTQ_FILE | $SAMTOOLS 'sort' -@ $THREADS > $BAM_FILE
$SAMTOOLS 'index' $BAM_FILE

$PILON --genome $FASTA_OUTPUT2 --nanopore $BAM_FILE --output ${SAMPLE}_pilon --outdir $OUT_DIR/${SAMPLE}_pilon --changes --chunksize 250000

# Start short-read polish: call polypolish pipeline
echo ""
echo -e "Step 6: Start short-read polish: polypolish pipeline"
PILON_FASTA=$OUT_DIR/${SAMPLE}_pilon/${SAMPLE}_pilon.fasta

$BWA index $PILON_FASTA
$BWA mem -t $THREADS -a $PILON_FASTA $PE1 > $OUT_DIR/${SAMPLE}_pilon/alignments_1.sam 
$BWA mem -t $THREADS -a $PILON_FASTA $PE2 > $OUT_DIR/${SAMPLE}_pilon/alignments_2.sam
$POLY_INSERT --in1 $OUT_DIR/${SAMPLE}_pilon/alignments_1.sam --in2 $OUT_DIR/${SAMPLE}_pilon/alignments_2.sam --out1 $OUT_DIR/${SAMPLE}_pilon/filtered_1.sam --out2 $OUT_DIR/${SAMPLE}_pilon/filtered_2.sam
$POLYPOLISH $PILON_FASTA $OUT_DIR/${SAMPLE}_pilon/filtered_1.sam $OUT_DIR/${SAMPLE}_pilon/filtered_2.sam > $OUT_DIR/${SAMPLE}_polypolish.fasta

echo ""
echo -e "Third rotation of assembly using random gene with dnaapler\n"
# Set the paths and names for input and output files for dnaapler2
FASTA_INPUT3=$OUT_DIR/${SAMPLE}_"polypolish.fasta"
FASTA_OUTPUT3=$OUT_DIR/${SAMPLE}_"dnaapler3.fasta"

# Call dnaapler mystery for third rotation of draft assembly
$SCRIPT_DIR/dnaapler_script.sh "$FASTA_INPUT3" "$FASTA_OUTPUT3" "$FUNCTION1" "$OUT_DIR"

echo ""
echo -e "Step 6: Start second short-read polish: Let's POLCA! - Derived from MaSuRCA-v4.1.0\n"

$POLCA -t $THREADS -a $FASTA_OUTPUT3 -r "$PE1 $PE2"

mkdir -p $OUT_DIR/${SAMPLE}_polca
mv $OUT_DIR/${SAMPLE}_dnaapler3_polca.fasta $OUT_DIR/${SAMPLE}_polca
mv $OUT_DIR/${SAMPLE}_polca/${SAMPLE}_dnaapler3_polca.fasta $OUT_DIR/${SAMPLE}_polca/${SAMPLE}_polca.fasta
mv $OUT_DIR/${SAMPLE}_dnaapler3_polca.report $OUT_DIR/${SAMPLE}_polca.report
mv $OUT_DIR/${SAMPLE}_polca.report/$OUT_DIR/logs

# Call dnaapler to reorient isolates one last time - run plasmid first given that plasmids will likely not have dnaA, thus will reorient on second call 

#Set the paths and names for input and output files for dnaapler3 and 4
FASTA_INPUT4=$OUT_DIR/${SAMPLE}_polca/${SAMPLE}_polca.fasta
FASTA_OUTPUT4=$OUT_DIR/${SAMPLE}_"dnaapler4.fasta"
FASTA_OUTPUT5=$OUT_DIR/${SAMPLE}_"consensus_flyest.fasta"
FUNCTION2="plasmid"
FUNCTION3="chromosome"

# Call dnaapler mystery for second rotation of draft assembly
echo ""
echo -e "Step 7: Final rotation of assembly using dnaA and repA genes with dnaapler\n"
$SCRIPT_DIR/dnaapler_script.sh "$FASTA_INPUT4" "$FASTA_OUTPUT4" "$FUNCTION2" "$OUT_DIR"
$SCRIPT_DIR/dnaapler_script.sh "$FASTA_OUTPUT4" "$FASTA_OUTPUT5" "$FUNCTION3" "$OUT_DIR"

echo ""
if [ "$EXECUTE_QUICK_QC" = true ]; then
    echo -e "Step 8: Executing quick_qc_script.sh to provide species, sequence type, contig count, contig length, and average coverage depth\n"
    # Call quick_qc_script.sh and pass the necessary parameters
    source $SCRIPT_DIR/quick_qc_script.sh $FASTA_OUTPUT5 $FASTQ_FILE $OUT_DIR $SAMPLE $THREADS
else
    echo -e "Skipping quick_qc_script.sh as --no-qc is specified.\n"
fi

# Move pipeline directories into a specified folder
mkdir -p $OUT_DIR/tmp
mkdir -p $OUT_DIR/${SAMPLE}_rotations
mv $OUT_DIR/*polypolish.fasta* $OUT_DIR/tmp
mv $OUT_DIR/*dnaapler* $OUT_DIR/${SAMPLE}_rotations
cp $OUT_DIR/${SAMPLE}_flye_assembly/assembly_info.txt $OUT_DIR/logs/${SAMPLE}_draft_assembly_info.txt
cp $OUT_DIR/${SAMPLE}_flye_assembly/flye.log $OUT_DIR/logs/${SAMPLE}_flye.log
mv $OUT_DIR/${SAMPLE}_polca $OUT_DIR/${SAMPLE}_rotations $OUT_DIR/${SAMPLE}_flye_assembly $OUT_DIR/${SAMPLE}_medaka $OUT_DIR/${SAMPLE}_berokka_clean $OUT_DIR/${SAMPLE}_pilon $OUT_DIR/tmp

# Remove any *bam files in pipeline directories to ensure minimal disk space constraint
echo -e "Remove large alignment files\n"
find $OUT_DIR/tmp -type f \( -name "*.bam" -o -name "*.bam.bai" -o -name "*sam" -o -name "*.sam*" \) -exec rm -r {} +

# Remove temporary files after pipeline is complete - will enable once at stable version of script
#if [ "$KEEP_TMP" = false]; then
#       rm -rf $OUT_DIR/tmp
#fi

# Record the end time
END_TIME=$(date +%s)

# Calculate the computation time in seconds
TOTAL_TIME=$((END_TIME - START_TIME))

# Convert total_time to days, hours, minutes, and seconds
DAYS=$((TOTAL_TIME / 86400))
TOTAL_TIME=$((TOTAL_TIME % 86400)) 
HOURS=$((TOTAL_TIME / 3600))
TOTAL_TIME=$((TOTAL_TIME % 3600)) 
MINUTES=$((TOTAL_TIME / 60))
SECONDS=$((TOTAL_TIME % 60))

# Print the result
echo "Finished"
echo "Computation time: $DAYS days, $HOURS hours, $MINUTES minutes, $SECONDS seconds"
echo "Pretty Flye-st for a genome! - Goodbye"
