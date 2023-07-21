#!/bin/bash

#Set the environment variables
export FASTA_FILE="$1"  # Path to your multi-fasta file
export OUTPUT_FILE="$2"  # Output file path
export FUNCTION="$3"	#dnaapler function
export OUTPUT_DIR="$4"  #output directory
output_directory="$OUTPUT_DIR/output"  # Output directory to store split files
dnaapler_dir="$OUTPUT_DIR/dnaapler"

# Create the output directory if it doesn't exist
mkdir -p "$output_directory"
mkdir -p "$dnaapler_dir"

echo "Separating multi-fasta file into individual files to be compliant with dnaapler"
# Read the multi-fasta file
while IFS= read -r line; do
  if [[ $line =~ ^\> ]]; then
    # Extract the sequence ID from the header line
    seq_id="${line#>}"

    # Create the directory with the appended number
    sub_directory="${output_directory}/${seq_id}"
    mkdir -p "$sub_directory"

    # Create the output file path
    output_file="${sub_directory}/${seq_id}.fasta"

    # Write the sequence header to the output file
    echo "$line" > "$output_file"
  else
    # Append the sequence to the output file
    echo "$line" >> "$output_file"
  fi
done < "$FASTA_FILE"

# Perform dnaapler on each file
count=1

for file in "$output_directory"/*/*.fasta; do
        # Perform dnaapler
        dnaapler_output_dir="$dnaapler_dir/SAMPLE_$count"
        dnaapler $FUNCTION -i "$file" -o "$dnaapler_output_dir" -p "SAMPLE_$count" -t 8

        #Check if fasta file is generated by dnaapler
        if find "$dnaapler_output_dir" -maxdepth 1 -type f -name "*.fasta" -print -quit | grep -q .; then
                echo "Performing dnaapler on $file"
        else
                echo "No fasta file generated by dnaapler for $file. Copying original input."
                cp "$file" "$dnaapler_output_dir/input.fasta"

        fi

        count=$((count + 1))
done

# Concatenate the individual files back into a single multi-fasta file
cat "$dnaapler_dir"/SAMPLE_*/*fasta > "$OUTPUT_FILE"

rm -rf $output_directory $dnaapler_dir