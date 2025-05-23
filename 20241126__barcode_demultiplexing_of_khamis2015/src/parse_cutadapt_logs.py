#!/usr/bin/env python3

import re
import glob
import os
import argparse
from pathlib import Path

def parse_cutadapt_log(log_file):
    """Parse a cutadapt log file and extract key metrics."""
    with open(log_file, 'r') as f:
        content = f.read()
    
    # Extract sample name from filename
    sample = Path(log_file).stem.replace('_cutadapt', '')
    
    # Extract key metrics using regex
    metrics = {
        'sample': sample,
        'total_pairs': int(re.search(r'Total read pairs processed:\s+([0-9,]+)', content).group(1).replace(',', '')),
        'pairs_with_adapter': int(re.search(r'Read 1 with adapter:\s+([0-9,]+)', content).group(1).replace(',', '')),
        'pairs_written': int(re.search(r'Pairs written \(passing filters\):\s+([0-9,]+)', content).group(1).replace(',', '')),
        'percent_written': float(re.search(r'Pairs written \(passing filters\):\s+[0-9,]+ \(([0-9.]+)%\)', content).group(1))
    }
    
    barcode_match = re.search(r'-g \^([ACGT]{6})', content)
    metrics['barcode'] = barcode_match.group(1) if barcode_match else "Unknown"
    
    return metrics

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Parse cutadapt log files and generate summary.')
    parser.add_argument('input_dir', help='Directory containing cutadapt log files')
    parser.add_argument('output_file', help='Path to output summary file')
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(args.output_file), exist_ok=True)
    
    # Find all cutadapt log files
    log_files = glob.glob(f"{args.input_dir}/*_cutadapt.log")
    
    if not log_files:
        print(f"No cutadapt log files found in {args.input_dir}!")
        return
    
    # Parse all log files
    results = []
    for log_file in log_files:
        try:
            metrics = parse_cutadapt_log(log_file)
            results.append(metrics)
        except Exception as e:
            print(f"Error parsing {log_file}: {str(e)}")
    
    # Sort results by sample name
    results.sort(key=lambda x: (x['sample'][0], int(x['sample'][1:]) if x['sample'][1:].isdigit() else float('inf')))
    
    # Write summary to file
    with open(args.output_file, 'w') as f:
        # Write header
        headers = ['Sample', 'Barcode', 'Total_Pairs', 'Pairs_With_Adapter', 'Pairs_Written', 'Percent_Written']
        f.write('\t'.join(headers) + '\n')
        
        # Write data
        for r in results:
            line = [
                r['sample'],
                r['barcode'],
                str(r['total_pairs']),
                str(r['pairs_with_adapter']),
                str(r['pairs_written']),
                f"{r['percent_written']:.1f}"
            ]
            f.write('\t'.join(line) + '\n')
    
    print(f"Summary written to: {args.output_file}")

if __name__ == "__main__":
    main()