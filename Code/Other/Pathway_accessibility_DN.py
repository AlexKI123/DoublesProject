import pandas as pd


# Dictionary mapping codon sequences to amino acids
# '_' indicates a stop codon
DNA_Codons = {
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
    "TAA": "Z", "TAG": "Z", "TGA": "Z"
}


def create_all_pathways(current_codon, alt_codon):
    return [
        current_codon[:i] + alt_codon[i] + current_codon[i+1:]
        for i in range(len(current_codon)) if current_codon[i] != alt_codon[i]
    ]


def get_fitness(codon, position, df, Wildtype):
    AAcurrent = DNA_Codons[Wildtype[position - 1]]
    AAmutant = DNA_Codons[codon]

    # Filter the DataFrame for the specific mutation
    mutation_df = df[df['Allele'] == "{}{}{}".format(AAcurrent, position, AAmutant)]

    # Check if there is a fitness value available
    if not mutation_df.empty:
        return float(mutation_df.iloc[0]['fitness'])
    else:
        return 1  # No AA change



def get_codons_for_amino_acid(aa):
    """Return codons corresponding to a given amino acid."""
    # Dictionary mapping amino acids to their codons
    codon_table = {
        'Z': ['TAA', 'TAG', 'TGA'],
        'A': ['GCT', 'GCC', 'GCA', 'GCG'],
        'C': ['TGT', 'TGC'],
        'D': ['GAT', 'GAC'],
        'E': ['GAA', 'GAG'],
        'F': ['TTT', 'TTC'],
        'G': ['GGT', 'GGC', 'GGA', 'GGG'],
        'H': ['CAT', 'CAC'],
        'I': ['ATT', 'ATC', 'ATA'],
        'K': ['AAA', 'AAG'],
        'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
        'M': ['ATG'],
        'N': ['AAT', 'AAC'],
        'P': ['CCT', 'CCC', 'CCA', 'CCG'],
        'Q': ['CAA', 'CAG'],
        'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
        'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
        'T': ['ACT', 'ACC', 'ACA', 'ACG'],
        'V': ['GTT', 'GTC', 'GTA', 'GTG'],
        'W': ['TGG'],
        'Y': ['TAT', 'TAC']
    }

    # Return codons for the provided amino acid
    return codon_table.get(aa, None)


def is_tandem_double_mutation(current_codon, alt_codons):
    """
    Check if any codons in alt_codons differ from current_codon by a tandem double nucleotide mutation.

    Args:
    - current_codon (str): The current codon.
    - alt_codons (list): List of alternative codons.

    Returns:
    - List of codons from alt_codons that differ by a tandem double nucleotide mutation from current_codon.
    """
    double_mutation_codons = []

    for alt_codon in alt_codons:

        # Identify positions where the codons differ
        diff_positions = [i for i, (a, b) in enumerate(zip(current_codon, alt_codon)) if a != b]

        # Check for tandem double mutations
        if len(diff_positions) == 2 and diff_positions[1] - diff_positions[0] == 1:
            double_mutation_codons.append(alt_codon)

    return double_mutation_codons


#data=pd.read_csv("/Users/alex/polybox/MNM/GithubDoubles/Data_P53/data_P53.txt", sep=",")
data=pd.read_csv("/Volumes/Dateien/ETHMAC/polybox_Copy/MNM/GithubDoubles/Data_P53/data_P53.txt", sep=",")


selected_tandem_data = data[data['mutation'] == "tandem"]
selected_tandem_data = selected_tandem_data[selected_tandem_data['GENIE_Mutation_Counts'] > 0]

#print(data['GENIE_Mutation_Counts'])
#print(data['fitness'])
print(selected_tandem_data)


# DNA sequence data
# P53
Wildtype = ['ATG', 'GAG', 'GAG', 'CCG', 'CAG', 'TCA', 'GAT', 'CCT', 'AGC', 'GTC', 'GAG', 'CCC', 'CCT', 'CTG', 'AGT', 'CAG', 'GAA', 'ACA', 'TTT', 'TCA', 'GAC', 'CTA', 'TGG', 'AAA', 'CTA', 'CTT', 'CCT', 'GAA', 'AAC', 'AAC', 'GTT', 'CTG', 'TCC', 'CCC', 'TTG', 'CCG', 'TCC', 'CAA', 'GCA', 'ATG', 'GAT', 'GAT', 'TTG', 'ATG', 'CTG', 'TCC', 'CCG', 'GAC', 'GAT', 'ATT', 'GAA', 'CAA', 'TGG', 'TTC', 'ACT', 'GAA', 'GAC', 'CCA', 'GGT', 'CCA', 'GAT', 'GAA', 'GCT', 'CCC', 'AGA', 'ATG', 'CCA', 'GAG', 'GCT', 'GCT', 'CCC', 'CGC', 'GTG', 'GCC', 'CCT', 'GCA', 'CCA', 'GCA', 'GCT', 'CCT', 'ACA', 'CCG', 'GCG', 'GCC', 'CCT', 'GCA', 'CCA', 'GCC', 'CCC', 'TCC', 'TGG', 'CCC', 'CTG', 'TCA', 'TCT', 'TCT', 'GTC', 'CCT', 'TCC', 'CAG', 'AAA', 'ACC', 'TAC', 'CAG', 'GGC', 'AGC', 'TAC', 'GGT', 'TTC', 'CGT', 'CTG', 'GGC', 'TTC', 'TTG', 'CAT', 'TCT', 'GGG', 'ACA', 'GCC', 'AAG', 'TCT', 'GTG', 'ACT', 'TGC', 'ACG', 'TAC', 'TCC', 'CCT', 'GCC', 'CTC', 'AAC', 'AAG', 'ATG', 'TTT', 'TGC', 'CAA', 'CTG', 'GCC', 'AAG', 'ACC', 'TGC', 'CCT', 'GTG', 'CAG', 'CTG', 'TGG', 'GTT', 'GAT', 'TCC', 'ACA', 'CCC', 'CCG', 'CCC', 'GGC', 'ACC', 'CGC', 'GTC', 'CGC', 'GCC', 'ATG', 'GCC', 'ATC', 'TAC', 'AAG', 'CAG', 'TCA', 'CAG', 'CAC', 'ATG', 'ACG', 'GAG', 'GTT', 'GTG', 'AGG', 'CGC', 'TGC', 'CCC', 'CAC', 'CAT', 'GAG', 'CGC', 'TGC', 'TCA', 'GAT', 'AGC', 'GAT', 'GGT', 'CTG', 'GCC', 'CCT', 'CCT', 'CAG', 'CAT', 'CTT', 'ATC', 'CGA', 'GTG', 'GAA', 'GGA', 'AAT', 'TTG', 'CGT', 'GTG', 'GAG', 'TAT', 'TTG', 'GAT', 'GAC', 'AGA', 'AAC', 'ACT', 'TTT', 'CGA', 'CAT', 'AGT', 'GTG', 'GTG', 'GTG', 'CCC', 'TAT', 'GAG', 'CCG', 'CCT', 'GAG', 'GTT', 'GGC', 'TCT', 'GAC', 'TGT', 'ACC', 'ACC', 'ATC', 'CAC', 'TAC', 'AAC', 'TAC', 'ATG', 'TGT', 'AAC', 'AGT', 'TCC', 'TGC', 'ATG', 'GGC', 'GGC', 'ATG', 'AAC', 'CGG', 'AGG', 'CCC', 'ATC', 'CTC', 'ACC', 'ATC', 'ATC', 'ACA', 'CTG', 'GAA', 'GAC', 'TCC', 'AGT', 'GGT', 'AAT', 'CTA', 'CTG', 'GGA', 'CGG', 'AAC', 'AGC', 'TTT', 'GAG', 'GTG', 'CGT', 'GTT', 'TGT', 'GCC', 'TGT', 'CCT', 'GGG', 'AGA', 'GAC', 'CGG', 'CGC', 'ACA', 'GAG', 'GAA', 'GAG', 'AAT', 'CTC', 'CGC', 'AAG', 'AAA', 'GGG', 'GAG', 'CCT', 'CAC', 'CAC', 'GAG', 'CTG', 'CCC', 'CCA', 'GGG', 'AGC', 'ACT', 'AAG', 'CGA', 'GCA', 'CTG', 'CCC', 'AAC', 'AAC', 'ACC', 'AGC', 'TCC', 'TCT', 'CCC', 'CAG', 'CCA', 'AAG', 'AAG', 'AAA', 'CCA', 'CTG', 'GAT', 'GGA', 'GAA', 'TAT', 'TTC', 'ACC', 'CTT', 'CAG', 'ATC', 'CGT', 'GGG', 'CGT', 'GAG', 'CGC', 'TTC', 'GAG', 'ATG', 'TTC', 'CGA', 'GAG', 'CTG', 'AAT', 'GAG', 'GCC', 'TTG', 'GAA', 'CTC', 'AAG', 'GAT', 'GCC', 'CAG', 'GCT', 'GGG', 'AAG', 'GAG', 'CCA', 'GGG', 'GGG', 'AGC', 'AGG', 'GCT', 'CAC', 'TCC', 'AGC', 'CAC', 'CTG', 'AAG', 'TCC', 'AAA', 'AAG', 'GGT', 'CAG', 'TCT', 'ACC', 'TCC', 'CGC', 'CAT', 'AAA', 'AAA', 'CTC', 'ATG', 'TTC', 'AAG', 'ACA', 'GAA', 'GGG', 'CCT', 'GAC', 'TCA', 'GAC', 'TAG']

position = 161
alt_codon = "F"

print("Wildtype codon: ", Wildtype[position-1])
#print("Wildtype AA: ", DNA_Codons[Wildtype[position-1]])
#print(get_codons_for_amino_acid("P"))

potential_tandem_mutations = is_tandem_double_mutation(Wildtype[position-1], get_codons_for_amino_acid(alt_codon))
print("potential_tandem_mutations", potential_tandem_mutations)
print("\n")
for i in range(len(potential_tandem_mutations)):
    Final_AA = potential_tandem_mutations[i]
    Final_AA_Fitness = get_fitness(Final_AA, position, data, Wildtype)
    print("Final_AA_Fitness", Final_AA_Fitness)
    intermediates = create_all_pathways(Wildtype[position-1], Final_AA)
    print("intermediates", intermediates)
    intermediate_1_fitness = get_fitness(intermediates[0], position, data, Wildtype)
    print("intermediate_1_fitness", intermediate_1_fitness)
    intermediate_2_fitness = get_fitness(intermediates[1], position, data, Wildtype)
    print("intermediate_2_fitness", intermediate_2_fitness)
    print("\n")

print(DNA_Codons['AGA'])
print(DNA_Codons['GAA'])
#accessible/mean fitness gain