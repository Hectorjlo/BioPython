# Importing the class of the previous work
from LópezOrdaz_HéctorJesús_Ejercicio2 import dna_sequence

###-------------------------### 
#   Ejercicio 3:
#   Usando el concepto de Herencia usando la clase gene, 
#   crea una llamada tRNA, y otra clase llamada RNA no codificante. 
#   Luego deriva de tRNA otra subclase llamada proteina
###-------------------------###


class tRNA(dna_sequence):
    has_cca_motif = False
    def check_cca(self):
        if len(self.sequence) >= 3:
            if self.sequence[-3:] == 'CCA':
                has_cca_motif = True
                result_analysis = "The tRNA has a CCA motif in the end 3'"
        else:
            result_analysis = 'Sequence too short to have a CCA motif'

        return result_analysis
        
        
 
        

        
        
       
    
if __name__ == '__main__':
    tRNA_sequence = tRNA('AUG-tRNA', 'GAUUGAGTACCG', 'normal expressed')
    print(f'El tRNA tiene la secuencia {tRNA_sequence.sequence} is attached to a {tRNA_sequence.name}')
    tRNA_sequence.check_cca()
    print()
    
    



    