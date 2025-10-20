#!/usr/bin/env bash

# Is meant to do the bwa indexing and the align at once


# Bariables
GENOME=""
PREFIX_INDEX="BWA"
INFILE_1=""
INFILE_2=""
OUTPUT=""
BATCH_MODE=false
REGEX_NAME=""
PARALLEL="1"


# Usage funtion
show_usage() {
    echo "Usage: $0 -i1 <input_file_1> -i2 <input_file_2> -o1 <output_file>"
    echo "   or: $0 -b -g <regex_pattern> -x <prefix_index> -ref <genome> -o <output_dir> [-p <parallel>]"
    echo ""
    echo "Required Arguments (single mode):"
    echo "  -ref,  --genome                  reference genome"
    echo "  -x,  --prefix                  prefix of the index (default: BWA)"
    echo "  -i1, --infile_1                paired file 1"
    echo "  -i2, --infile_2                paired file 2"
    echo "  -o,  --output                  output file" 
    echo "Example:"
    echo ""
    echo "Batch Mode:"
    echo "  -b,   --batch                  enable batch processing mode"
    echo "  -g,   --regex_name             regex pattern to find paired files (e.g., '*_1.fastq.gz')"
    echo "  -ref, --genome                 reference genome"
    echo "  -x,   --prefix                 prefix of the index (default: BWA)"
    echo "  -o,   --output                 output directory for batch results"
    echo "  -p,   --parallel               parallel execution (default: 1)"
    echo ""
    echo "Example (single):"
    echo "  $0 -ref reference.fasta.gz -x BWA -i1 trimmed_1.fastq.gz -i2 trimmed_2.fastq.gz -o aligned.sam"
    echo "Example (batch):"
    echo "  $0 -b -g 'trimmed*_1.fastq.gz' -ref reference.fasta.gz -x BWA -o aligned_results -p 4"
    exit 1
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -b|--batch)
            BATCH_MODE=true
            shift 1
            ;;
        -g|--regex_name)
            REGEX_NAME="$2"
            shift 2
            ;;
        -ref|--genome)
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
        -p|--parallel)
            PARALLEL="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            show_usage
            ;;
    esac
done


# Funtion to aling a single pair of files
align_single() {
    local genome="$1"
    local prefix="$2"
    local infile_1="$3"
    local infile_2="$4"
    local output="$5"

    echo ""
    echo "[INFO] Aligning: $(basename "$infile_1")"
    echo "  Reference: $genome"
    echo "  Input R1:  $infile1"
    echo "  Input R2:  $infile2"
    echo "  Output:    $output"
    echo ""

    # Lauch bwa mem
    bwa mem -o "$output" "$prefix" "$infile1" "$infile2"

    #verify if completed
    if [ $? -eq 0 ]; then
        echo "[DONE] Completed: $(basename "$output")"
    else
        echo "[ERROR] Failed: $(basename "$output")"
        return 1
    fi

}

process_batch() {
    local infile_1="$1"

    # Get the name and the dir
    namebase=$(basename "$infile_1")
    dir=$(dirname "$infile_1")

    #Get paired 2
    infile_2="${dir}/${namebase/_1.fastq.gz/_2.fastq.gz}"
    # Check if exists
    if [ ! -f "$infile_2" ]; then
        echo "[ERROR] Paired file not found: $infile_2"
        return 1
    fi

    output_file="${OUTPUT}/aligned_${namebase/_1.fastq.gz/.sam}"

    echo ""
    align_single "$GENOME" "$PREFIX_INDEX" "$infile_1" "$infile_2" "$output_file"


}

# condition of batch mode and verify the completness of the variables
if [ "$BATCH_MODE" = true ]; then
    # Batch mode validation
    if [[ -z "$GENOME" || -z "$REGEX_NAME" || -z "$OUTPUT" ]]; then
        echo "Error: Batch mode requires -ref, -g, and -o arguments"
        show_usage
    fi
    

    # Create output directory if needed
    mkdir -p "$OUTPUT"

    # Index genome once (if not already indexed)
    echo "[INFO] Indexing genome: $GENOME"
    bwa index -p "$PREFIX_INDEX" "$GENOME"
    echo "" 

    # Export for parallel
    export -f align_single
    export -f process_batch
    export GENOME PREFIX_INDEX OUTPUT

    # Process files in parallel
    find . -path "${REGEX_NAME}" | parallel -j "$PARALLEL" process_batch {}

else
    # Single mode if not:
    # Check the variables
    if [[ -z "$GENOME" || -z "$INFILE_1" || -z "$INFILE_2" || -z "$OUTPUT" ]]; then
        echo "Error: All arguments are required for single mode"
        show_usage
    fi

    # Create output if needed
    output_dir=$(dirname "$OUTPUT")
    mkdir -p "$output_dir"

    # Index and align
    echo "[INFO] Indexing genome: $GENOME"
    bwa index -p "$PREFIX_INDEX" "$GENOME"
    echo ""

    # funtion
    align_single "$GENOME" "$PREFIX_INDEX" "$INFILE_1" "$INFILE_2" "$OUTPUT"

fi


echo ""
echo ""
echo "[DONE] Thank you for using align_bwa!"