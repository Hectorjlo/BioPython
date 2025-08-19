# Importing the class of the previous work
from LópezOrdaz_HéctorJesús_Ejercicio2 import gene_sequence

###-------------------------### 
#   Ejercicio 3:
#   Usando el concepto de Herencia usando la clase gene, 
#   crea una llamada tRNA, y otra clase llamada RNA no codificante. 
#   Luego deriva de tRNA otra subclase llamada proteina
###-------------------------###


class tRNA(gene_sequence):

    # Separate by 3 pairs given a reading frame default with 1
    def split_by_codons(self, reading_frame = 1):
        codons = 1 #[for codon in gene_sequence.sequence] 

if __name__ == '__main__':
    dna = gene_sequence('lacY', 'AAATGCATCAGCTAGCTAGCTAGCATCGATCGAT', 'overexpressed')

    print(dna.sequence)


    