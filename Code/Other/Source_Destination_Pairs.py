# Define the standard genetic code
genetic_code = {
    # 'M' - START, '_' - STOP
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TGT": "C", "TGC": "C",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "TTT": "F", "TTC": "F",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAT": "H", "CAC": "H",
    "ATA": "I", "ATT": "I", "ATC": "I",
    "AAA": "K", "AAG": "K",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATG": "M",
    "AAT": "N", "AAC": "N",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TGG": "W",
    "TAT": "Y", "TAC": "Y",
    "TAA": "_", "TAG": "_", "TGA": "_"
}

# Initialize an empty set to store the amino acid changes
amino_acid_changes = set()

# List of nucleotides
nucleotides = ['A', 'C', 'G', 'T']

# Generate all possible codons
codons = [n1 + n2 + n3 for n1 in nucleotides for n2 in nucleotides for n3 in nucleotides]

# Iterate through each codon
for codon in codons:
    # Skip stop codons
    if genetic_code[codon] == '_':
        continue

    # Generate all possible single nucleotide mutations
    for i in range(3):
        for mutated_nucleotide in nucleotides:
            if codon[i] != mutated_nucleotide:
                mutated_codon = codon[:i] + mutated_nucleotide + codon[i+1:]
                # Skip mutations that result in stop codons
                if genetic_code[mutated_codon] == '_':
                    continue
                # Check if the mutation changes the amino acid
                if genetic_code[codon] != genetic_code[mutated_codon]:
                    amino_acid_changes.add((genetic_code[codon], genetic_code[mutated_codon]))

# Print the size of the set
print(f"Size of the set of codon to amino acid changes: {len(amino_acid_changes)}")
