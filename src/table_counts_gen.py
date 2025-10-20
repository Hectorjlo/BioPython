#!/usr/bin/env python3
"""
Generate a counts table from multiple featureCounts output files.
Combines multiple samples into a single table with gene names and read counts.
"""

import pandas as pd
import argparse
import re
from pathlib import Path
import sys
from typing import List, Optional, Tuple


def parse_gene_info(info_string: str) -> Tuple[str, str]:
    """
    Extract gene ID and name from the annotation column.
    
    Args:
        info_string: String like 'ID=FBgn0267594;Name=CR45932;Alias=CG40583;...'
    
    Returns:
        Tuple of (gene_id, gene_name)
    """
    gene_id = ""
    gene_name = ""
    
    # Extract ID
    id_match = re.search(r'ID=([^;]+)', info_string)
    if id_match:
        gene_id = id_match.group(1)
    
    # Extract Name
    name_match = re.search(r'Name=([^;]+)', info_string)
    if name_match:
        gene_name = name_match.group(1)
    
    return gene_id, gene_name


def read_featurecounts_file(filepath: str, sample_name: Optional[str] = None) -> pd.DataFrame:
    """
    Read a featureCounts output file and extract relevant columns.
    
    Args:
        filepath: Path to the featureCounts file
        sample_name: Optional name for the sample (defaults to filename)
    
    Returns:
        DataFrame with columns: gene_id, gene_name, reads, covered_bases, gene_length, coverage_pct
    """
    if sample_name is None:
        sample_name = Path(filepath).stem
    
    data = []
    
    with open(filepath, 'r') as f:
        for line in f:
            # Skip comment lines
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            
            # Expecting at least 13 columns
            if len(fields) < 13:
                continue
            
            # Parse gene info from column 9 (index 8)
            gene_id, gene_name = parse_gene_info(fields[8])
            
            # Extract counts and coverage info
            try:
                reads = int(fields[9])
                covered_bases = int(fields[10])
                gene_length = int(fields[11])
                coverage_pct = float(fields[12])
                
                data.append({
                    'gene_id': gene_id,
                    'gene_name': gene_name,
                    'sample': sample_name,
                    'reads': reads,
                    'covered_bases': covered_bases,
                    'gene_length': gene_length,
                    'coverage_pct': coverage_pct
                })
            except ValueError:
                # Skip lines with non-numeric values
                continue
    
    return pd.DataFrame(data)


def create_counts_table(input_files: List[str], sample_names: Optional[List[str]] = None, output_file: Optional[str] = None) -> pd.DataFrame:
    """
    Create a combined counts table from multiple featureCounts files.
    
    Args:
        input_files: List of paths to featureCounts files
        sample_names: Optional list of sample names (same order as input_files)
        output_file: Optional output file path (if None, prints to stdout)
    
    Returns:
        DataFrame with combined counts
    """
    all_data = []
    
    # Read all files
    for i, filepath in enumerate(input_files):
        sample_name = sample_names[i] if sample_names else None
        print(f"[INFO] Reading {filepath}...", file=sys.stderr)
        df = read_featurecounts_file(filepath, sample_name)
        all_data.append(df)
    
    # Combine all data
    combined_df = pd.concat(all_data, ignore_index=True)
    
    # Create pivot table: rows = genes, columns = samples, values = reads
    counts_table = combined_df.pivot_table(
        index=['gene_id', 'gene_name'],
        columns='sample',
        values='reads',
        fill_value=0
    )
    
    # Reset index to make gene_id and gene_name regular columns
    counts_table = counts_table.reset_index()
    
    # Save or print
    if output_file:
        counts_table.to_csv(output_file, sep='\t', index=False)
        print(f"[DONE] Counts table saved to {output_file}", file=sys.stderr)
    else:
        print(counts_table.to_csv(sep='\t', index=False))
    
    return counts_table


def main() -> None:
    """
    Main entry point for the script. Parses arguments and generates the counts table.
    """
    parser = argparse.ArgumentParser(
        description='Generate counts table from featureCounts output files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Single file
  python table_counts_gen.py -i counts1.txt -o counts_table.tsv
  
  # Multiple files with custom names
  python table_counts_gen.py -i counts1.txt counts2.txt counts3.txt \\
                             -n sample1 sample2 sample3 \\
                             -o counts_table.tsv
  
  # Using wildcards
  python table_counts_gen.py -i counts/*.txt -o counts_table.tsv
        '''
    )
    
    parser.add_argument('-i', '--input', nargs='+', required=True,
                        help='Input featureCounts file(s)')
    parser.add_argument('-n', '--names', nargs='+',
                        help='Sample names (same order as input files)')
    parser.add_argument('-o', '--output',
                        help='Output file (TSV format). If not provided, prints to stdout')
    
    args = parser.parse_args()
    
    # Validate input
    if args.names and len(args.names) != len(args.input):
        parser.error("Number of sample names must match number of input files")
    
    # Check if files exist
    for filepath in args.input:
        if not Path(filepath).exists():
            print(f"[ERROR] File not found: {filepath}", file=sys.stderr)
            sys.exit(1)
    
    # Create counts table
    create_counts_table(args.input, args.names, args.output)


if __name__ == "__main__":
    main()