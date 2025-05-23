#!/bin/bash

# alorenzetti 20241126

# description ####
# this script will take the fastq files
# and demultiplex them based on barcodes

# setup ####
# Check if required programs are available
if ! command -v python3 &> /dev/null; then
    echo "Error: python3 not found"
    exit 1
fi

# cutadapt path
cutadapt_path="/opt/homebrew/Caskroom/miniforge/base/envs/condaenv__amellifera/bin/cutadapt"

# Check if cutadapt is available
if ! command -v ${cutadapt_path} &> /dev/null; then
    echo "Error: cutadapt not found at ${cutadapt_path}"
    exit 1
fi

# set number of cores to use
threads=8

# getting started ####
function demultiplex_fastq() {
    local input_r1=$1
    local input_r2=$2
    local barcode=$3
    local sample=$4
    # Update output paths to use output_dir
    local output_r1="${output_dir}/${sample}_1.fastq.gz"
    local output_r2="${output_dir}/${sample}_2.fastq.gz"
    
    # According to the cutadapt documentation
    # the --pair-filter=any will be automatically replaced
    # by --pair-filter=both for this specific case (-g and --discard-untrimmed)
    ${cutadapt_path} \
        -j ${threads} \
        -g "^${barcode}NNNNNNNNTATAGGG" \
        --discard-untrimmed \
        --pair-filter=any \
        -o "${output_r1}" \
        -p "${output_r2}" \
        "${input_r1}" \
        "${input_r2}" \
        > "${output_dir}/${sample}_cutadapt.log" \
        2> "${output_dir}/${sample}_cutadapt.err"
    
    # Check cutadapt exit status
    if [ $? -ne 0 ]; then
        echo "Error: Cutadapt failed for sample ${sample}"
        exit 1
    fi
    
    echo "Cutadapt finished for sample: $sample"
    echo "Check ${sample}_cutadapt.log and ${sample}_cutadapt.err for details"
}

# Barcodes
barcodes1_keys=(ACAGAT ATCGTG CACGAT CACTGA CTGACG GAGTGA GTATAC TCGAGC)
barcodes1_values=(N2 N4 N29 N31 N32 N33 N41 N42)

barcodes2_keys=(ACAGAT ATCGTG CACGAT CACTGA CTGACG GAGTGA GTATAC TCGAGC)
barcodes2_values=(F14 F27 F28 F41 F48 F49 F50 F54)

# Validate barcode arrays have matching lengths
if [ ${#barcodes1_keys[@]} -ne ${#barcodes1_values[@]} ]; then
    echo "Error: Mismatch in barcodes1 arrays length"
    exit 1
fi

if [ ${#barcodes2_keys[@]} -ne ${#barcodes2_values[@]} ]; then
    echo "Error: Mismatch in barcodes2 arrays length"
    exit 1
fi

# Create output directory if it doesn't exist
output_dir="data/processed"
if [ ! -d "$output_dir" ]; then
    mkdir -p "$output_dir"
fi

# making a minimal example #
# Create test subset of fastq files using pigz for faster compression
# unpigz -p ${threads} -c data/raw/SRR1582004_1.fastq.gz | head -n 100000 | pigz -p ${threads} > data/raw/SRR1582004_1.test.fastq.gz
# unpigz -p ${threads} -c data/raw/SRR1582004_2.fastq.gz | head -n 100000 | pigz -p ${threads} > data/raw/SRR1582004_2.test.fastq.gz

# Update the loop to use array indices
echo "Starting demultiplexing process..."

# running a minimal example
# for i in "${!barcodes1_keys[@]}"; do
#     barcode="${barcodes1_keys[$i]}"
#     sample="${barcodes1_values[$i]}"
#     echo "Processing barcode: $barcode"
#     echo "Sample name: $sample"
    
#     # Check if input files exist
#     if [ ! -f "data/raw/SRR1582004_1.test.fastq.gz" ] || [ ! -f "data/raw/SRR1582004_2.test.fastq.gz" ]; then
#         echo "Error: Input files not found"
#         exit 1
#     fi
    
#     demultiplex_fastq "data/raw/SRR1582004_1.test.fastq.gz" "data/raw/SRR1582004_2.test.fastq.gz" "$barcode" "$sample"
# done

echo "Available barcodes: ${barcodes1_keys[@]}"
for i in "${!barcodes1_keys[@]}"; do
    barcode="${barcodes1_keys[$i]}"
    sample="${barcodes1_values[$i]}"
    echo "Processing barcode: $barcode" 
    echo "Sample name: $sample"
    
    # Check if input files exist
    if [ ! -f "data/raw/SRR1582004_1.fastq.gz" ] || [ ! -f "data/raw/SRR1582004_2.fastq.gz" ]; then
        echo "Error: Input files not found"
        exit 1
    fi
    
    demultiplex_fastq "data/raw/SRR1582004_1.fastq.gz" "data/raw/SRR1582004_2.fastq.gz" "$barcode" "$sample"
done

echo "Available barcodes: ${barcodes2_keys[@]}"
for i in "${!barcodes2_keys[@]}"; do
    barcode="${barcodes2_keys[$i]}"
    sample="${barcodes2_values[$i]}"
    echo "Processing barcode: $barcode"
    echo "Sample name: $sample"
    
    # Check if input files exist
    if [ ! -f "data/raw/SRR1582203_1.fastq.gz" ] || [ ! -f "data/raw/SRR1582203_2.fastq.gz" ]; then
        echo "Error: Input files not found"
        exit 1
    fi
    
    demultiplex_fastq "data/raw/SRR1582203_1.fastq.gz" "data/raw/SRR1582203_2.fastq.gz" "$barcode" "$sample"
done

echo "Demultiplexing process completed"

# Parse cutadapt logs and generate summary table
echo "Generating summary of cutadapt results..."
python3 src/parse_cutadapt_logs.py "data/processed/" "data/processed/cutadapt_demultiplexing_summary.tsv"
echo "Summary table generated at data/processed/cutadapt_demultiplexing_summary.tsv"
