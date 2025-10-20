#!/usr/bin/env bash

# The following script is meant to use cut adapt to clean various .fastq files


# Paths for the files
REGEX_NAME=""
QUALITY_CUTOFF=""
MIN_LENGTH=""
OUTPUT_PREFIX="trimmed"
CORES="1"
PARALLEL="1"


# Usage funtion
show_usage() {
    echo "Usage: $0 -i1 <input_file_1> -i2 <input_file_2> -o1 <output_file_1> -o2 <output_file_2>"
    echo ""
    echo "Required Arguments:"
    echo "  -g,  --regex_name              name in regex to search recursivelly"
    echo "  -q,  --quality-cutoff          cutadapt -q"
    echo "  -m,  --minimum-length          cutadapt -m"
    echo "  -o,  --output_prefix           cutadapt -o (default: trimmed)"
    echo "  -j,  --cores                   cutadapt -j (default: 1), actual use would be the cores times parallel value"
    echo "  -p,  --parallel                parallel execution (default: 1)"
    echo "Example:"
    echo "  $0 -g 'sample*.fastq.gz' -q 20 -m 25 -p trimmed -j 10"
    exit 1
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -g|--regex_name)
            REGEX_NAME="$2"
            shift 2
            ;;
        -q|--quality-cutoff)
            QUALITY_CUTOFF="$2"
            shift 2
            ;;
        -m|--minimum-length)
            MIN_LENGTH="$2"
            shift 2
            ;;
        -o|--output_prefix)
            OUTPUT_PREFIX="$2"
            shift 2
            ;;
        -j|--cores)
            CORES="$2"
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

# Check if all required arguments are provided
if [[ -z "$REGEX_NAME" || -z "$QUALITY_CUTOFF" || -z "$MIN_LENGTH" || -z "$OUTPUT_PREFIX" || -z "$CORES" ]]; then
    echo "Error: All arguments are required"
    show_usage
fi


# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_PREFIX}"

# Funtion to run cut adapt
cut_adapters() {
    infile_1="$1"

    # Get base name and the dir
    namebase=$(basename "$infile_1")
    dir=$(dirname "$infile_1")

    # Using the sustitutions get paired 2
    infile_2="${dir}/${namebase/_1.fastq.gz/_2.fastq.gz}"

    # Cechk if that file exists
    if [ ! -f "$infile_2" ]; then
        echo "[ERROR] Paired file not found: $infile_2"
        return 1
    fi

    # Output files and the log
    outfile_1="${OUTPUT_PREFIX}/${OUTPUT_PREFIX}_${namebase}"
    outfile_2="${OUTPUT_PREFIX}/${OUTPUT_PREFIX}_${namebase/_1.fastq.gz/_2.fastq.gz}"
    log_file="${OUTPUT_PREFIX}/${OUTPUT_PREFIX}_${namebase}.log"


    echo ""
    echo "[INFO] Processing: $namebase"
    echo "  Input R1:  $infile_1"
    echo "  Input R2:  $infile_2"
    echo "  Output R1: $outfile_1"
    echo "  Output R2: $outfile_2"
    echo ""

    # Run cutadapt
    cutadapt -q "$QUALITY_CUTOFF" -m "$MIN_LENGTH" -j "$CORES" --pair-filter=any -o "$outfile_1" -p "$outfile_2" "$infile_1" "$infile_2" > "$log_file" 2>&1
    
    if [ $? -eq 0 ]; then
        echo "[DONE] Completed: $namebase"
    else
        echo "[ERROR] Failed: $namebase (check log: $log_file)"
        return 1
    fi
}


# Export function and variables for parallel
export -f cut_adapters
export QUALITY_CUTOFF MIN_LENGTH OUTPUT_PREFIX CORES


find . -name "${REGEX_NAME}" | parallel -j "$PARALLEL" cut_adapters {}
echo ""
echo ""
echo "[DONE] Thank you for using run_cutadapter_pairend!"