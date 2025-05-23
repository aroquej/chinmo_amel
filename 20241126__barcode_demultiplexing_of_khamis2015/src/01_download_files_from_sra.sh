#!/bin/bash

# alorenzetti 20241126

# description ####
# this script will download the files from the SRA database
# using links obtained from SRA-Explorer for the given accession number
# PRJNA261549, Khamis2015, https://www.nature.com/articles/srep11136

# setup ####
# the requirements for this script to work are:
# sra-tools (prefetch, fastq-dump)

# the directory where the files will be downloaded
raw_data_dir="data/raw"
if [ ! -d $raw_data_dir ]; then
    mkdir -p $raw_data_dir
fi

# getting started ####
# SRR1582004 is nurses
# SRR1582203 is foragers

# First prefetch the SRA files
for accession in SRR1582004 SRR1582203; do
    echo "Prefetching $accession..."
    if ! prefetch -O $raw_data_dir $accession; then
        echo "Error prefetching $accession"
        exit 1
    fi
done

# Then convert to fastq
for accession in SRR1582004 SRR1582203; do
    echo "Converting $accession to fastq..."
    if ! fastq-dump --split-files --gzip --outdir $raw_data_dir "${raw_data_dir}/${accession}/${accession}.sra"; then
        echo "Error converting $accession to fastq"
        exit 1
    fi
done

echo "All files downloaded and converted successfully"