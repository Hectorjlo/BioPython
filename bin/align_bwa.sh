#!/usr/bin/env bash

# Is meant to do the bwa indexing and the align at once


# Bariables
GENOME=""
PREFIX_INDEX="BWA"
INFILE_1=""
INFILE_2=""
OUTPUT=""


# Usage funtion
show_usage() {
    echo "Usage: $0 -i1 <input_file_1> -i2 <input_file_2> -o1 <output_file_1> -o2 <output_file_2>"
    echo ""
    echo "Required Arguments:"
    echo "  -g,  --genome                  reference genome"
    echo "  -x,  --prefix                  prefix of the index"
    echo "  -i1, --infile_1                paired file 1"
    echo "  -i2, --infile_2                paired file 2"
    echo "  -o,  --output                  output file" 
    echo "Example:"
    echo "  $0 -g 'reference.fasta.gz' -x BWA -i1 trimmed_1.fastq.gz -i2 trimmed_2.fastq.gz -p trimmed"
    exit 1
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -g|--genome)
            GENOME="$2"
            shift 2
            ;;
        -x|--prefix)
            PREFIX_INDEX="$2"
            shift 2
            ;;
        -i1|--infile_1)
            INFILE_1="$2"
            shift 2
            ;;
        -i2|--infile_2)
            INFILE_2="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            show_usage
            ;;
    esac
done

# Check if all required arguments are provided
if [[ -z "$GENOME" || -z "$INFILE_1" || -z "$INFILE_2" || -z "$OUTPUT" ]]; then
    echo "Error: All arguments are required"
    show_usage
fi


# Create output directory if needed
output_dir=$(dirname "$OUTPUT")
mkdir -p "$output_dir"  

# Run the two at once

bwa index -p "$PREFIX_INDEX" "$GENOME"

bwa mem -o "$OUTPUT" "$PREFIX_INDEX" "$INFILE_1" "$INFILE_2"


echo ""
echo ""
echo "[DONE] Thank you for using align_bwa!"