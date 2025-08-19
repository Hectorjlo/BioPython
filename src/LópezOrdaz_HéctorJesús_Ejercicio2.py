class gene_sequence():

    def __init__(self, name, sequence, expression):
        self.name = name
        self.expression = expression
        self.sequence = sequence.upper().replace(" ", "")


    def lenght(self):
        """Gives the lenght of the sequence"""
        return len(self.sequence)

    def gc_content(self):
        """GC content of the sequence"""
        if not self.sequence:
            return 0.0
        g = self.sequence.count('G')
        c = self.sequence.count('C')
        return round((g + c) * 100.0 / len(self.sequence), 2)

    
# Probe 

if __name__ == '__main__':

    dna_gene = gene_sequence("lacY", "AAATGCATCAGCTAGCTAGCTAGCATCGATCGAT", "overexpressed")
    print(f'The lenght of the gene ({dna_gene.name}) is of {dna_gene.lenght()} bases and is {dna_gene.expression}')
    print(f'It has {dna_gene.gc_content()}% of GC content')