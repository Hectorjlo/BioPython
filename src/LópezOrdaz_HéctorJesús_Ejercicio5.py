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

def seq_obj_creator(sequence: str) -> Seq:
    """
    Transforms a str sequence and a str sense in q Seq object and return it with the sense str

    Arguments:
        sequence (str): The sequence in str format, DNA or RNA
        
    Returns:
        obj_sequence (Seq): Seq object of Bio
    """
    # Hence in the problem the sense is not given, I'll be using like a 5' --> 3' 
    sequence_input = sequence.strip()

    # Creating the object seq()
    obj_sequence = Seq(sequence_input)

    return obj_sequence


def trancribe_seq_obj(obj_sequence: Seq, sequence_type: str = 'DNA') -> Seq:
    
    # Give info to the user 
    sequence_legth = len(obj_sequence)
    if sequence_legth > 1:
        print(f"[INFO] The sequence has a legth of: {sequence_legth} {'nucleotides' if sequence_legth > 1 else 'nucleotide'}")
    elif sequence_legth <= 0:
        print(f"[ERROR] The sequence legth is a invalid value")

    # It needed to transform DNA to RNA

    if sequence_type == 'DNA':
        # Using .transcribe()
        transcribed_sequence = obj_sequence.transcribe()
    # If is already RNA
    elif sequence_type == 'RNA':
        # No need to transcribe
        transcribed_sequence = obj_sequence

    return transcribed_sequence

def get_orfs(transcribed_sequence: Seq, start_codon: str|list = 'AUG', is_complement: bool = False) -> list:

    # Initalize a where orfs will save
    orfs = []
    # Using a sequence it will check all the 3 ORFs

    # Initalize the start and stop codons
    start_codon = 'AUG'
    stop_codons = ['UAA', 'UAG', 'UGA']

    # Look for the start codon in each position of the sequence to check all 3 frames
    for reading_frame in range(3):
        frame_label = reading_frame + 1 if not is_complement else (reading_frame + 4)
        print(f'[INFO] Proccessing frame: {frame_label}')
        # Check in steps of 3 by 3
        for i in range(reading_frame, len(transcribed_sequence), 3):
            check_codon = str(transcribed_sequence[i: i+3])
            # Check if a match is made
            if check_codon == start_codon:
                orf_start = i
                # Check from previous match to the end of the sequence or stop codon
                for j in range(i, len(transcribed_sequence), 3):
                    # Check the following triplets 
                    codon_proccesed = str(transcribed_sequence[j: j+3])
                    # If a stop codon is found or the codons are insufficient to check a triplet
                    if codon_proccesed in stop_codons or len(codon_proccesed) < 3:
                        # Extract the ORF
                        orf_sequence = transcribed_sequence[orf_start: j]
                        if len(orf_sequence) >= 3:
                            orfs.append((orf_sequence, f'length:{len(orf_sequence)}', f'Reading Frame:{frame_label}'))
                        break
    
    return orfs

    



def find_all_reading_frames(transcribed_sequence: Seq) -> Seq:
    # There are 6 reading frames, 3 in the forward strand but other 3 on the reverse strand
    # Get the orfs
    orfs_first_sense = get_orfs(transcribed_sequence)

    # Get the reverse complement 
    reversed_complement = transcribed_sequence.reverse_complement()
      
    orfs_second_sense = get_orfs(reversed_complement, is_complement=True)

    final_orfs = orfs_first_sense + orfs_second_sense
    
    return final_orfs

def translate_orfs(orfs: list):
    proteins = []

    for orf in orfs:
        protein = orf[0].translate(to_stop=True)
        proteins.append([protein, orf[1], orf[2]])
    
    return proteins

    
   





        
def main() -> None:
    obj_sequence = seq_obj_creator(sequence='AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG')
    transcribed_sequence = trancribe_seq_obj(obj_sequence=obj_sequence, sequence_type='DNA')
    orfs = find_all_reading_frames(transcribed_sequence=transcribed_sequence)
    proteins = translate_orfs(orfs=orfs)

    print("\n[DONE] Protein sequences found:")
    for i, protein_info in enumerate(proteins, 1):
        print(f"\nProtein {i}:")
        print(f"    Sequence: {protein_info[0]}")
        print(f"    {protein_info[1]}, {protein_info[2]}")
    
    # Find longest protein
    if proteins:
        longest = max(proteins, key=lambda x: len(x[0]))
        print("\n[DONE] Longest protein sequence:")
        print(f"    Sequence: {longest[0]}")
        print(f"    {longest[1]}, {longest[2]}")

    
if __name__ == '__main__':
    main()
