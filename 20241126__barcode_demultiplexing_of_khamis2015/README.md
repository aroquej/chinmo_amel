# Khamis et al. 2015 Honeybee RNA-Seq demultiplexing

## Overview
This pipeline processes RNA-Seq data from honeybee (Apis mellifera) nurses and foragers, performing barcode-based demultiplexing of CAGE samples. The data comes from the study by Khamis et al. (2015) published in Scientific Reports (https://www.nature.com/articles/srep11136).

## Data Source
- **SRA Project**: PRJNA261549
- **Samples**:
  - SRR1582004 (Nurses)
  - SRR1582203 (Foragers)

## Barcodes
The samples were multiplexed using template switching oligos (TSOs) having 6-bp barcodes on their 5' end:

### Nurses (SRR1582004)
| Sample | Barcode |
|--------|---------|
| N2     | ACAGAT  |
| N4     | ATCGTG  |
| N29    | CACGAT  |
| N31    | CACTGA  |
| N32    | CTGACG  |
| N33    | GAGTGA  |
| N41    | GTATAC  |
| N42    | TCGAGC  |

### Foragers (SRR1582203)
| Sample | Barcode |
|--------|---------|
| F14    | ACAGAT  |
| F27    | ATCGTG  |
| F28    | CACGAT  |
| F41    | CACTGA  |
| F48    | CTGACG  |
| F49    | GAGTGA  |
| F50    | GTATAC  |
| F54    | TCGAGC  |


See Table S8 in Khamis et al. 2015 for more details.

## Pipeline Structure

### 1. Data Download
The script `src/01_download_files_from_sra.sh` downloads the raw sequencing data from SRA using `prefetch` and `fastq-dump` tools.

### 2. Demultiplexing
The script `src/02_demultiplex_fastq_files.sh` processes the multiplexed FASTQ files using cutadapt to separate individual samples based on their barcodes. The script:
- Processes paired-end reads
- Identifies and trims barcodes
- Generates separate FASTQ files for each sample
- Creates detailed log files for each demultiplexing

### 3. Summary Generation
A Python script (`src/parse_cutadapt_logs.py`) analyzes the cutadapt log files to generate a summary table containing:
- Sample names
- Barcodes
- Total read pairs
- Pairs with adapters
- Pairs written
- Percentage of reads retained

## Requirements
- SRA Tools (fastq-dump, prefetch)
- cutadapt
- pigz
- Python 3.x

## Usage
1. Ensure all requirements are installed
2. Run the driver script:
   ```bash
   bash src/run_pipeline.sh
   ```

## Results
The pipeline generates:
1. Individual FASTQ files for each sample in `data/processed/`
2. Detailed cutadapt logs for each demultiplexing
3. A summary table (`cutadapt_deconvolution_summary.tsv`) showing the success rate of barcode identification and sample separation

## License
This project is licensed under the MIT License - see the LICENSE file for details.