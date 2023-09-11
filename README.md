# Flyest

## Nanopore Only Pipeline - Flyer | Nanopore + Illumina Hybrid Consensus Polishing Pipeline - Flyest - README

## Objective
This pipeline was created for the purpose of polishing bacterial draft genome assemblies from the Flye assembler using either Oxford Nanopore Technologies long-reads only (Flyer) or a hybrid approach that combines ONT long-read and Illumina short-read sequencing data. The pipeline integrates multiple tools and performs iterative rotations of the assembly to improve the quality of the consensus genome. This script is inspired by looking at best practices as of 2023-07-17 for both long-read only and hybrid assembly approaches with a lot of inspiration from the work of Ryan Wick. 

## Author

[William Shropshire](https://twitter.com/The_Real_Shrops)

## Scope
The pipeline covers the following steps:
1. Initial draft assembly using Flye with long-read input (Nanopore).
2. Evaluation and filtering of contigs based on minimum coverage depth and minimum contig length.
3. First rotation of the draft assembly using dnaapler.
4. Medaka polish using a specified Medaka model for error correction.
5. Berokka trimming of overlaps if present in conjunction with a modified clean script to remove identical contigs.
6. Second rotation of the assembly using dnaapler.
7. Final long-read polish using Pilon.
8. Optional: Short-read polish using polypolish pipeline followed by assembly roation and a second short-read polish using an adapted polca.sh script if Flyest is enabled.
9. Final rotation of the assembly using dnaapler for proper contig orientation using dnaA and repA genes.
10. Quick QC to determine putative bacterial species, sequence type, contig number, contig length, and mean coverage depth of consensus assembly. 

## Installation and Dependencies
The pipeline requires the following programs to be available in the system's PATH:
- Flye
- Dnaapler 
- Medaka
- Pilon
- Berokka
- BWA
- Mummer4
- Polypolish (Flyest dependency)
- MLST

The pipeline utilizes a custom script `clean.py` adapted from `circlator-v1.5.5` to remove highly similar contigs based on Nucmer identity and lengths. The `polca_mod.sh` short-read polishing script is adapted from MaSuRCA-v4.1.0 to work within the flyest environment and output files in a user created output directory within the pipeline. I would highly recommend creating a conda environment and using the conda pack tool with the **flyest_v0.1.tar.gz** package available in this GitHub. If installing manually, a list of all dependencies is included in the **flyest_v0.1_conda_env.yml** file. The `ufasta` binary is in the **binaries** directory. 

Also note that the Conda package only includes a limited number of Medaka models and that institutional firewalls may impede from downloading models if running on an HPC. I recommend downloading their models from their data directory and placing model files in the appropriate directory. 

## Function and Usage
The script takes various command-line arguments to control its behavior:

```
Usage: flyest_pipeline.sh [options]

Options:
  -i, --input       File containing long-read fastq reads [Can be gzipped] -- REQUIRED.
  -1, --pe1         Paired-end forward short-reads -- REQUIRED.
  -2, --pe2         Paired-end reverse short-reads -- REQUIRED.
  --nano-corr       Long-read input quality [Default=--nano-hq for ONT high-quality reads; Only use --nano-corr for pre Guppy-v5.x basecalled reads].
  -n, --min         Minimum coverage depth per sequence as a fraction of the overall draft assembly mean coverage depth [Default=0.2; i.e. 20%]. Use a value of '0.001' for no minimum coverage depth.
  -l, --len         Minimum contig length per sequence [Default=1000].
  --meta            Meta option for Flye [Default is not set] - Good for uneven/low coverage assemblies.
  -m, --mod         Medaka model [Default=r1041_e82_400bps_sup_v4.2.0].
  --nid             Nucmer min_id parameter for removing highly similar contigs based on percent minimum nucleotide identity shared [Default = 95].
  --nlen            Nucmer min_length parameter for removing highly similar contigs based on percent minimum contig length shared [Default = 90].
  --no-qc           Disable quick QC script [Default is on].
  -o, --outdir      Output directory -- REQUIRED.
  -t, --threads     Number of threads to use [Default=1].
  -s, --sample      Sample prefix [DEFAULT=SAMPLE].
  -h, --help        Show this help text.
```

The only difference for Flyer is that paired-end short-read inputs are not an option: 

```
Usage: flyest_pipeline.sh [options]

Options:
  -i, --input       File containing long-read fastq reads [Can be gzipped] -- REQUIRED.
  --nano-corr       Long-read input quality [Default=--nano-hq for ONT high-quality reads; Only use --nano-corr for pre Guppy-v5.x basecalled reads].
  -n, --min         Minimum coverage depth per sequence as a fraction of the overall draft assembly mean coverage depth [Default=0.2; i.e. 20%]. Use a value of '0.001' for no minimum coverage depth.
  -l, --len         Minimum contig length per sequence [Default=1000].
  --meta            Meta option for Flye [Default is not set] - Good for uneven/low coverage assemblies.
  -m, --mod         Medaka model [Default=r1041_e82_400bps_sup_v4.2.0].
  --nid             Nucmer min_id parameter for removing highly similar contigs based on percent minimum nucleotide identity shared [Default = 95].
  --nlen            Nucmer min_length parameter for removing highly similar contigs based on percent minimum contig length shared [Default = 90].
  --no-qc           Disable quick QC script [Default is on].
  -o, --outdir      Output directory -- REQUIRED.
  -t, --threads     Number of threads to use [Default=1].
  -s, --sample      Sample prefix [DEFAULT=SAMPLE].
  -h, --help        Show this help text.
```

## Usage Example
```bash
flyest_pipeline.sh -i long_reads.fastq.gz -1 pe_forward.fastq.gz -2 pe_reverse.fastq.gz -o output_dir --meta -t 16 -s SAMPLE_NAME
```

## Note
The script assumes that necessary software and reference files are available or configured in the system, including a Medaka model and necessary genome assembly files.
