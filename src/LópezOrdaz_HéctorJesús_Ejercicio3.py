# Importing the class of the previous work
from LópezOrdaz_HéctorJesús_Ejercicio2 import dna_sequence

###-------------------------### 
#   Ejercicio 3:
#   Usando el concepto de Herencia usando la clase gene, 
#   crea una llamada tRNA, y otra clase llamada RNA no codificante. 
#   Luego deriva de tRNA otra subclase llamada proteina
###-------------------------###

## Hola Mario :) :
# Puede que no fuera el mejor código, pero me limité a no cambiar tanto el primer ejercicio
# se me hizo complicado el hecho de contruir los objetos al paso, y no tener la idea clara desde el inicio
# principalmente en los argumentos del primer objeto y que los demás toman, por eso tanto el uso de los returns
# y la complejidad de que biologicamente era complicado hilar todo entre si, pero si me quedo claro el uso de la herencia
# pero estuvo bien



class tRNA(dna_sequence):
    
    def check_cca(self):
        if len(self.sequence) >= 3:
            if self.sequence[-3:] == 'CCA':
                result_analysis = "The tRNA has a CCA motif in the end 3'"
                return result_analysis
            else:
                result_analysis = "The tRNA has not a CCA motif in the end 3'"
                return result_analysis
        else:
            result_analysis = 'Sequence too short to have a CCA motif'
            return result_analysis
    
    def anticodon_check(self, anticodon):
        if anticodon.upper() in self.sequence:
            anticodon_motif = "Found"
            return anticodon_motif
        else:
            anticodon_motif = "Not found"
            return anticodon_motif
    # I use length_ (I notice that I don't spell it right but ._. )
    def lenght_(self):
        return self.lenght()


class ncRNA(dna_sequence):

    def size_check(self):
        lenght = self.lenght()
        if lenght >= 200:
            rna_type = "lncRNA"
            return rna_type
        else:
            rna_type = 'sncRNA'
            return rna_type

class protein(tRNA):
    #Using overriding
    def __init__(self, name, sequence, expression, dna_input=""):
        tRNA.__init__(self, name, sequence, expression)
        self.dna_input = dna_input
        
    def calculate_molecular_weight(self):
        aa_weights = {
            'A': 89.1, 'R': 174.2, 'N': 132.1, 'D': 133.1, 'C': 121.2,
            'E': 147.1, 'Q': 146.2, 'G': 75.1, 'H': 155.2, 'I': 131.2,
            'L': 131.2, 'K': 146.2, 'M': 149.2, 'F': 165.2, 'P': 115.1,
            'S': 105.1, 'T': 119.1, 'W': 204.2, 'Y': 181.2, 'V': 117.1,
        }
        
        molecular_weight = sum(aa_weights.get(aminoacid, 0.0) for aminoacid in self.sequence) 
        
        return molecular_weight
    
    def protein_lenght(self):
        aminoacid_lenght = self.lenght()
        return (aminoacid_lenght)
    
    def nucleotide_lenght(self):
        if self.dna_input != "":
            nucleotide_lenght = len(self.dna_input)
        else:
            nucleotide_lenght = "(!Not provided)"

        return nucleotide_lenght
    


        
       
    
if __name__ == '__main__':
    tRNA_sequence = tRNA('AUG-tRNA', 'GAUUGAGTACCG', 'normal expressed')
    print(f'The sequence of the tRNA is: {tRNA_sequence.sequence}')
    result_analysis = tRNA_sequence.check_cca()
    print(f"After the analysis of the tRNA 3' end it's showed that ' {result_analysis} '")
    anticodon_motif = tRNA_sequence.anticodon_check("GAU")
    print(f'The anticodon provided was {anticodon_motif}')
    rna_type = ncRNA("let-7", 'UAGCUAGUCGCUAGCUGAUCGAUGCUAGCUAGUCGAUCGUAGCUAGCUAGUCAUGCUAGCUAGCUGAUCGUAGCUAGUCGAUGCUAGCUGAUCGAUCGUAGCUAGUCGAUCGUAGUGAUGCUAGUCGAUCGUAGUUUUCGACGACAGCGACAGGCAGCUAGCAACACGAGCAGCGAGCGAUCGAUGCUAGCUAGUCGAUGC', 'normal-expressed')
    print(f'The name of the ncRNA is: {rna_type.name}')
    size_check = rna_type.size_check()
    print(f"By the lenght the ncRNA is a: '{size_check}'")
    protein_sequence = protein('LacY', 'MYYLKNTNFWMFGLFFFFYFFIMGAYFPFFPIWLHDINHISKSDTGIIFAAISLFSLLFQPLFGLLSDKLGLRKYLLWIITGMLVMFAPFFIFIFGPLLQYNILVGSIVGGIYLGFCFNAGAPAVEAFIEKVSRRSNFEFGRARMFGCVGWALCASIVGIMFTINNQFVFWLGSGCALILAVLLFFAKTDAPSSATVANAVGANHSAFSLKLALELFRQPKLWFLSLYVIGVSCTYDVFDQQFANFFTSFFATGEQGTRVFGYVTTMGELLNASIMFFAPLIINRIGGKNALLLAGTIMSVRIIGSSFATSALEVVILKTLHMFEVPFLLVGCFKYITSQFEVRFSATIYLVCFCFFKQLAMIFMSVLAGNMYESIGFQGAYLVLGLVALGFTLISVFTLSGPGPLSLLRRQVNEVA', 'Normal-expressed', 'ATGACCATGATTACGGATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGGTCGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCG')
    molecular_weight = protein_sequence.calculate_molecular_weight()
    print(f"The molecular weight of the protein '{protein_sequence.name}' is of {round(molecular_weight, 4)} Da")
    rna_lenght = tRNA_sequence.lenght_()
    print(f"The length of the tRNA is of: {rna_lenght} nucleotides")
    protein_lenght = protein_sequence.protein_lenght()
    protein_nucleotides_amount = protein_sequence.nucleotide_lenght()
    print(f"The total amount of aminoacids is of: {protein_lenght}, (If given) The total nucleotides were of: {protein_nucleotides_amount}")





    