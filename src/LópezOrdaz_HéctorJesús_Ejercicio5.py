"""
Ejercicio 5: 
    Obtener cadena protéica de cualquiera de sus ORFs
    Un ORF inicia con un codón inicial y termina, ya sea con un codón final o al final de la cadena. 

    Input a utilizar:

    AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG
    Ejercicio avanzado (opcional):
    Elegir la cadena protéica de mayor longitud (No olvidar probar todos los ORFs)

"""

# Importing Biopython 
from Bio.Seq import Seq

def seq_obj_creator(sequence: str) -> tuple[Seq, Seq]:
    """
    Create Seq objects for both forward and reverse complement strands.
    
    This function takes a DNA sequence string and generates two Seq objects:
    one for the complement of the forward strand and one for the reverse complement.
    Both are oriented 5' to 3'.
    
    Args:
        sequence (str): The input DNA sequence as a string (assumed 5' to 3' orientation)
        
    Returns:
        tuple[Seq, Seq]: A tuple containing:
            - obj_sequence_forward (Seq): Complement of the input sequence (5' to 3')
            - obj_sequence_reversed (Seq): Reverse complement of the input sequence (5' to 3')

    """
    # Remove any leading/trailing whitespace from input
    sequence_input = sequence.strip()

    # Creating the Seq objects for both strands
    obj_sequence_forward = Seq(sequence_input).complement()            # Complement: 5' to 3'
    obj_sequence_reversed = Seq(sequence_input).reverse_complement()  # Reverse complement: 5' to 3'

    return obj_sequence_forward, obj_sequence_reversed


def trancribe_seq_obj(obj_sequence_forward: Seq, obj_sequence_reversed: Seq) -> tuple[Seq, Seq]:
    """
    Transcribe DNA sequences to mRNA for both forward and reversed strands.
    
    This function converts DNA Seq objects to their corresponding mRNA sequences
    by replacing thymine (T) with uracil (U). It also validates the sequence length
    and provides informative output.
    
    Args:
        obj_sequence_forward (Seq): The forward DNA sequence object
        obj_sequence_reversed (Seq): The reversed DNA sequence object
        
    Returns:
        tuple[Seq, Seq]: A tuple containing:
            - mrna_forward (Seq): Transcribed mRNA from forward strand
            - mrna_reversed (Seq): Transcribed mRNA from reversed strand

    
    """
    # Validate sequence length and provide user feedback
    sequence_legth = len(obj_sequence_forward)
    if sequence_legth > 1:
        print(f"[INFO] The sequence has a legth of: {sequence_legth} {'nucleotides' if sequence_legth > 1 else 'nucleotide'}")
    elif sequence_legth <= 0:
        print(f"[ERROR] The sequence legth is a invalid value")

    # Transcribe DNA to mRNA (T -> U)
    mrna_forward = obj_sequence_forward.transcribe()
    mrna_reversed = obj_sequence_reversed.transcribe()

    return mrna_forward, mrna_reversed


def get_orfs(transcribed_sequence: Seq, start_codon: str = 'AUG', frame_offset: int = 0) -> list:
    """
    Find all Open Reading Frames (ORFs) in a given RNA sequence.
    
    An ORF is defined as a sequence that:
    1. Starts with a start codon (default: AUG)
    2. Ends with a stop codon (UAA, UAG, or UGA) or at the end of the sequence
    3. Has a length of at least 3 nucleotides
    
    This function analyzes all three reading frames of the input sequence.
    
    Args:
        transcribed_sequence (Seq): The mRNA sequence to analyze
        start_codon (str, optional): The start codon to search for. Defaults to 'AUG'
        frame_offset (int, optional): Offset for frame numbering. Use 0 for forward 
            strand (frames 1-3) and 3 for reverse strand (frames 4-6). Defaults to 0
    
    Returns:
        list: List of tuples, where each tuple contains:
            - orf_sequence (Seq): The ORF nucleotide sequence
            - length_info (str): String describing the ORF length (e.g., 'length:45')
            - frame_info (str): String describing the reading frame (e.g., 'Reading Frame:1')
    """
    # Initialize list where orfs will be saved
    orfs = []
    
    # Define stop codons for translation termination
    stop_codons = ['UAA', 'UAG', 'UGA']

    # Analyze all three reading frames
    for reading_frame in range(3):
        # Calculate frame label (1-3 for forward, 4-6 for reverse)
        frame_label = reading_frame + 1 + frame_offset
        
        print(f'[INFO] Processing frame: {frame_label}')
        start_codons_found = 0
        
        # Scan the sequence in steps of 3 nucleotides (codon by codon)
        for i in range(reading_frame, len(transcribed_sequence), 3):
            check_codon = str(transcribed_sequence[i: i+3])
            
            # When a start codon is found, look for the end of the ORF
            if check_codon == start_codon:
                start_codons_found += 1
                orf_start = i
                orf_found = False
                
                # Search from start codon to end of sequence or stop codon
                for j in range(orf_start, len(transcribed_sequence), 3):
                    codon_processed = str(transcribed_sequence[j: j+3])
                    
                    # Case 1: Stop codon found
                    if codon_processed in stop_codons:
                        orf_sequence = transcribed_sequence[orf_start: j]
                        if len(orf_sequence) >= 3:
                            orfs.append((orf_sequence, f'length:{len(orf_sequence)}', f'Reading Frame:{frame_label}'))
                        orf_found = True
                        break
                    # Case 2: Incomplete codon at end of sequence
                    elif len(codon_processed) < 3:
                        orf_sequence = transcribed_sequence[orf_start: j]
                        if len(orf_sequence) >= 3:
                            orfs.append((orf_sequence, f'length:{len(orf_sequence)}', f'Reading Frame:{frame_label}'))
                        orf_found = True
                        break
                
                # Case 3: ORF extends to the end of the sequence
                if not orf_found:
                    orf_sequence = transcribed_sequence[orf_start:]
                    if len(orf_sequence) >= 3:
                        orfs.append((orf_sequence, f'length:{len(orf_sequence)}', f'Reading Frame:{frame_label}'))
                        
    
    return orfs


def find_all_orfs(mrna_forward: Seq, mrna_reversed: Seq) -> list:
    """
    Find all ORFs in both forward and reverse complement RNA sequences.
    
    This function analyzes all six reading frames (three in the forward strand and
    three in the reverse complement strand) to identify all possible ORFs.
    
    Args:
        mrna_forward (Seq): Forward mRNA sequence (frames 1-3)
        mrna_reversed (Seq): Reverse complement mRNA sequence (frames 4-6)
    
    Returns:
        list: Combined list of all ORFs from both strands. Each element is a tuple
            containing (orf_sequence, length_info, frame_info)
    
    """
    # Analyze forward strand (reading frames 1-3)
    print("\n[INFO] Analyzing forward strand...")
    orfs_forward = get_orfs(mrna_forward, frame_offset=0)
    
    # Analyze reverse complement strand (reading frames 4-6)
    print("\n[INFO] Analyzing reversed strand...")
    orfs_reversed = get_orfs(mrna_reversed, frame_offset=3)
    
    # Combine results from both strands
    all_orfs = orfs_forward + orfs_reversed

  
    return all_orfs


def translate_orfs(orfs: list) -> list:
    """
    Translate ORF nucleotide sequences to their corresponding protein sequences.
    
    This function takes a list of ORFs and translates each one to its amino acid
    sequence using the standard genetic code. Translation stops at the first
    stop codon encountered.
    
    Args:
        orfs (list): List of tuples containing (orf_sequence, length_info, frame_info)
            where orf_sequence is a Seq object representing the nucleotide sequence
    
    Returns:
        list: List of lists, where each inner list contains:
            - protein (Seq): The translated protein sequence
            - length_info (str): Length information from the original ORF
            - frame_info (str): Reading frame information from the original ORF
    """
    proteins = []

    # Translate each ORF to protein sequence
    for orf in orfs:
        # Translate to stop codon (to_stop=True stops at first stop codon)
        protein = orf[0].translate(to_stop=True)
        proteins.append([protein, orf[1], orf[2]])
    
    return proteins

    
   





        
def main() -> None:
    obj_sequence_forward, obj_sequence_reversed = seq_obj_creator(sequence='AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG')
    mrna_forward, mrna_reversed = trancribe_seq_obj(obj_sequence_forward=obj_sequence_forward, obj_sequence_reversed=obj_sequence_reversed)
    orfs = find_all_orfs(mrna_forward=mrna_forward, mrna_reversed=mrna_reversed)
    
    # Step 4: Translate ORFs to protein sequences
    proteins = translate_orfs(orfs=orfs)

    # Step 5: Display all protein sequences found
    print("\n[DONE] Protein sequences found:")
    for i, protein_info in enumerate(proteins, 1):
        print(f"\nProtein {i}:")
        print(f"    Sequence: {protein_info[0]}")
        print(f"    Length: {int(protein_info[1].split(':')[1])/3}, {protein_info[2]}")
    
    # Step 6: Find and display the longest protein (advanced feature)
    if proteins:
        longest = max(proteins, key=lambda x: len(x[0]))
        print("\n[DONE] Longest protein sequence:")
        print(f"    Sequence: {longest[0]}")
        print(f"    Length: {int(longest[1].split(':')[1])/3}, {longest[2]}")

    
if __name__ == '__main__':
    main()
