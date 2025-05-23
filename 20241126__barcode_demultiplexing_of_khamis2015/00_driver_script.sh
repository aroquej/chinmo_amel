#!/bin/bash

# alorenzetti 20241126

# description ####
# this script is a driver script for the project
# it will source the other scripts and execute them in the correct order

# setup ####
# creating data folder
data_dir="data"
if [ ! -d $data_dir ]; then
    mkdir $data_dir
fi

# creating config folder
config_dir="config"
if [ ! -d $config_dir ]; then
    mkdir $config_dir
fi

# writing requirements file
requirements_file="requirements.txt"
if [ ! -f $requirements_file ]; then
    touch $requirements_file
fi

# writing readme file
readme_file="README.md"
if [ ! -f $readme_file ]; then
    touch $readme_file
fi

# adding requirements to the file
echo "fastq-dump" >> $requirements_file
echo "prefetch" >> $requirements_file
echo "cutadapt" >> $requirements_file
echo "pigz" >> $requirements_file

# getting started ####
bash src/01_download_files_from_sra.sh
bash src/02_demultiplex_fastq_files.sh
