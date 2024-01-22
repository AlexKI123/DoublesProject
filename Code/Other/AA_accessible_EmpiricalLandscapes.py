# Define DNA nucleotides as both a list and a set
DNA_Nucleotides = ['A', 'C', 'G', 'T']
DNA_Nucleotides_set = {'A', 'C', 'G', 'T'}

# Dictionary for DNA reverse complement
DNA_ReverseComplement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

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
    "TAA": "_", "TAG": "_", "TGA": "_"
}


# Function to compute the number of different characters between two strings
def differing_characters(s1, s2):
    return sum(1 for a, b in zip(s1, s2) if a != b)


# Function to convert a nucleotide sequence to an amino acid sequence
def AASequence(ntSeq):
    Codons = int(len(ntSeq)/3)
    AA = []
    for c in range(Codons):
        codonSeq = ''.join(ntSeq[0+3*c:3+3*c])
        AA.append(DNA_Codons[codonSeq])
    return ''.join(AA)



# DNA sequence data
# P53
#Wildtype = ['ATG', 'GAG', 'GAG', 'CCG', 'CAG', 'TCA', 'GAT', 'CCT', 'AGC', 'GTC', 'GAG', 'CCC', 'CCT', 'CTG', 'AGT', 'CAG', 'GAA', 'ACA', 'TTT', 'TCA', 'GAC', 'CTA', 'TGG', 'AAA', 'CTA', 'CTT', 'CCT', 'GAA', 'AAC', 'AAC', 'GTT', 'CTG', 'TCC', 'CCC', 'TTG', 'CCG', 'TCC', 'CAA', 'GCA', 'ATG', 'GAT', 'GAT', 'TTG', 'ATG', 'CTG', 'TCC', 'CCG', 'GAC', 'GAT', 'ATT', 'GAA', 'CAA', 'TGG', 'TTC', 'ACT', 'GAA', 'GAC', 'CCA', 'GGT', 'CCA', 'GAT', 'GAA', 'GCT', 'CCC', 'AGA', 'ATG', 'CCA', 'GAG', 'GCT', 'GCT', 'CCC', 'CGC', 'GTG', 'GCC', 'CCT', 'GCA', 'CCA', 'GCA', 'GCT', 'CCT', 'ACA', 'CCG', 'GCG', 'GCC', 'CCT', 'GCA', 'CCA', 'GCC', 'CCC', 'TCC', 'TGG', 'CCC', 'CTG', 'TCA', 'TCT', 'TCT', 'GTC', 'CCT', 'TCC', 'CAG', 'AAA', 'ACC', 'TAC', 'CAG', 'GGC', 'AGC', 'TAC', 'GGT', 'TTC', 'CGT', 'CTG', 'GGC', 'TTC', 'TTG', 'CAT', 'TCT', 'GGG', 'ACA', 'GCC', 'AAG', 'TCT', 'GTG', 'ACT', 'TGC', 'ACG', 'TAC', 'TCC', 'CCT', 'GCC', 'CTC', 'AAC', 'AAG', 'ATG', 'TTT', 'TGC', 'CAA', 'CTG', 'GCC', 'AAG', 'ACC', 'TGC', 'CCT', 'GTG', 'CAG', 'CTG', 'TGG', 'GTT', 'GAT', 'TCC', 'ACA', 'CCC', 'CCG', 'CCC', 'GGC', 'ACC', 'CGC', 'GTC', 'CGC', 'GCC', 'ATG', 'GCC', 'ATC', 'TAC', 'AAG', 'CAG', 'TCA', 'CAG', 'CAC', 'ATG', 'ACG', 'GAG', 'GTT', 'GTG', 'AGG', 'CGC', 'TGC', 'CCC', 'CAC', 'CAT', 'GAG', 'CGC', 'TGC', 'TCA', 'GAT', 'AGC', 'GAT', 'GGT', 'CTG', 'GCC', 'CCT', 'CCT', 'CAG', 'CAT', 'CTT', 'ATC', 'CGA', 'GTG', 'GAA', 'GGA', 'AAT', 'TTG', 'CGT', 'GTG', 'GAG', 'TAT', 'TTG', 'GAT', 'GAC', 'AGA', 'AAC', 'ACT', 'TTT', 'CGA', 'CAT', 'AGT', 'GTG', 'GTG', 'GTG', 'CCC', 'TAT', 'GAG', 'CCG', 'CCT', 'GAG', 'GTT', 'GGC', 'TCT', 'GAC', 'TGT', 'ACC', 'ACC', 'ATC', 'CAC', 'TAC', 'AAC', 'TAC', 'ATG', 'TGT', 'AAC', 'AGT', 'TCC', 'TGC', 'ATG', 'GGC', 'GGC', 'ATG', 'AAC', 'CGG', 'AGG', 'CCC', 'ATC', 'CTC', 'ACC', 'ATC', 'ATC', 'ACA', 'CTG', 'GAA', 'GAC', 'TCC', 'AGT', 'GGT', 'AAT', 'CTA', 'CTG', 'GGA', 'CGG', 'AAC', 'AGC', 'TTT', 'GAG', 'GTG', 'CGT', 'GTT', 'TGT', 'GCC', 'TGT', 'CCT', 'GGG', 'AGA', 'GAC', 'CGG', 'CGC', 'ACA', 'GAG', 'GAA', 'GAG', 'AAT', 'CTC', 'CGC', 'AAG', 'AAA', 'GGG', 'GAG', 'CCT', 'CAC', 'CAC', 'GAG', 'CTG', 'CCC', 'CCA', 'GGG', 'AGC', 'ACT', 'AAG', 'CGA', 'GCA', 'CTG', 'CCC', 'AAC', 'AAC', 'ACC', 'AGC', 'TCC', 'TCT', 'CCC', 'CAG', 'CCA', 'AAG', 'AAG', 'AAA', 'CCA', 'CTG', 'GAT', 'GGA', 'GAA', 'TAT', 'TTC', 'ACC', 'CTT', 'CAG', 'ATC', 'CGT', 'GGG', 'CGT', 'GAG', 'CGC', 'TTC', 'GAG', 'ATG', 'TTC', 'CGA', 'GAG', 'CTG', 'AAT', 'GAG', 'GCC', 'TTG', 'GAA', 'CTC', 'AAG', 'GAT', 'GCC', 'CAG', 'GCT', 'GGG', 'AAG', 'GAG', 'CCA', 'GGG', 'GGG', 'AGC', 'AGG', 'GCT', 'CAC', 'TCC', 'AGC', 'CAC', 'CTG', 'AAG', 'TCC', 'AAA', 'AAG', 'GGT', 'CAG', 'TCT', 'ACC', 'TCC', 'CGC', 'CAT', 'AAA', 'AAA', 'CTC', 'ATG', 'TTC', 'AAG', 'ACA', 'GAA', 'GGG', 'CCT', 'GAC', 'TCA', 'GAC', 'TAG']

# TEM 1
Wildtype = ['ATG', 'AGT', 'ATT', 'CAA', 'CAT', 'TTC', 'CGT', 'GTC', 'GCC', 'CTT', 'ATT', 'CCC', 'TTT', 'TTT', 'GCG', 'GCA', 'TTT', 'TGC', 'CTT', 'CCT', 'GTT', 'TTT', 'GCT', 'CAC', 'CCA', 'GAA', 'ACG', 'CTG', 'GTG', 'AAA', 'GTA', 'AAA', 'GAT', 'GCT', 'GAA', 'GAT', 'CAG', 'TTG', 'GGT', 'GCA', 'CGA', 'GTG', 'GGT', 'TAC', 'ATC', 'GAA', 'CTG', 'GAT', 'CTC', 'AAC', 'AGC', 'GGT', 'AAG', 'ATC', 'CTT', 'GAG', 'AGT', 'TTT', 'CGC', 'CCC', 'GAA', 'GAA', 'CGT', 'TTT', 'CCA', 'ATG', 'ATG', 'AGC', 'ACT', 'TTT', 'AAA', 'GTT', 'CTG', 'CTA', 'TGT', 'GGC', 'GCG', 'GTA', 'TTA', 'TCC', 'CGT', 'GTT', 'GAC', 'GCC', 'GGG', 'CAA', 'GAG', 'CAA', 'CTC', 'GGT', 'CGC', 'CGC', 'ATA', 'CAC', 'TAT', 'TCT', 'CAG', 'AAT', 'GAC', 'TTG', 'GTT', 'GAG', 'TAC', 'TCA', 'CCA', 'GTC', 'ACA', 'GAA', 'AAG', 'CAT', 'CTT', 'ACG', 'GAT', 'GGC', 'ATG', 'ACA', 'GTA', 'AGA', 'GAA', 'TTA', 'TGC', 'AGT', 'GCT', 'GCC', 'ATA', 'ACC', 'ATG', 'AGT', 'GAT', 'AAC', 'ACT', 'GCG', 'GCC', 'AAC', 'TTA', 'CTT', 'CTG', 'ACA', 'ACG', 'ATC', 'GGA', 'GGA', 'CCG', 'AAG', 'GAG', 'CTA', 'ACC', 'GCT', 'TTT', 'TTG', 'CAC', 'AAC', 'ATG', 'GGG', 'GAT', 'CAT', 'GTA', 'ACT', 'CGC', 'CTT', 'GAT', 'CGT', 'TGG', 'GAA', 'CCG', 'GAG', 'CTG', 'AAT', 'GAA', 'GCC', 'ATA', 'CCA', 'AAC', 'GAC', 'GAG', 'CGT', 'GAC', 'ACC', 'ACG', 'ATG', 'CCT', 'GCA', 'GCA', 'ATG', 'GCA', 'ACA', 'ACG', 'TTG', 'CGC', 'AAA', 'CTA', 'TTA', 'ACT', 'GGC', 'GAA', 'CTA', 'CTT', 'ACT', 'CTA', 'GCT', 'TCC', 'CGG', 'CAA', 'CAA', 'TTA', 'ATA', 'GAC', 'TGG', 'ATG', 'GAG', 'GCG', 'GAT', 'AAA', 'GTT', 'GCA', 'GGA', 'CCA', 'CTT', 'CTG', 'CGC', 'TCG', 'GCC', 'CTT', 'CCG', 'GCT', 'GGC', 'TGG', 'TTT', 'ATT', 'GCT', 'GAT', 'AAA', 'TCT', 'GGA', 'GCC', 'GGT', 'GAG', 'CGT', 'GGG', 'TCT', 'CGC', 'GGT', 'ATC', 'ATT', 'GCA', 'GCA', 'CTG', 'GGG', 'CCA', 'GAT', 'GGT', 'AAG', 'CCC', 'TCC', 'CGT', 'ATC', 'GTA', 'GTT', 'ATC', 'TAC', 'ACG', 'ACG', 'GGG', 'AGT', 'CAG', 'GCA', 'ACT', 'ATG', 'GAT', 'GAA', 'CGA', 'AAT', 'AGA', 'CAG', 'ATC', 'GCT', 'GAG', 'ATA', 'GGT', 'GCC', 'TCA', 'CTG', 'ATT', 'AAG', 'CAT', 'TGG']


Codons = len(Wildtype)
Loci = Codons*3
my_list = ''.join(Wildtype)
# Convert the string to a tuple of characters
wild_type = tuple(my_list)  # Original DNA sequence
AA = AASequence(wild_type)  # Amino acid sequence of the wild type


# Counters for various mutation types and their effects
counterSN = 0  # counter for AA accessibility through SN
counterDN = 0  # counter for AA accessibility through DN
counterDNMore = 0  # counter for AA accessibility through DN - SN
counterSNsyn = 0  # counter for synonymous SN mutation
counterDNsyn = 0  # counter for synonymous DN mutation

counterMutationsSN = 0

NewAASetSN = {AA}
# Iterate over each nucleotide position (locus) and consider single nucleotide mutations (SN)
for locus in range(Loci):  # iterating through all loci
    substi = list(DNA_Nucleotides_set - {wild_type[locus]})  # all substitutions
    copy_item = list(wild_type[:])
    for s in range(len(substi)):  # iterating through all substitutions
        counterMutationsSN +=1
        copy_item[locus] = substi[s]
        NewAA = AASequence(copy_item)  # mutant AA
        #if locus >= 391:
        #    print(NewAA[-1])
        NewAASetSN.add(NewAA)
        if AA == NewAA:  # Check for synonymous mutations
            counterSNsyn += 1
counterSN = len(NewAASetSN) - 1  # Count unique amino acid sequences due to SN mutations

print("\nAccessible Amino acids through SN mutations", counterSN)

counterMutationsDN = 0
counterMutationsDN_twoAAmutations = 0

# Generate all possible double nucleotide (DN) mutation combinations
dn_mutations = [[i,i+1] for i in range(Codons*3-1)]   # all DN combinations
NewAASetDN = {AA}  # wild type AA
NewTwoAASetDN = {AA}
NewOneAASetDN = {AA}

# Iterate over each DN mutation combination
for dn_comb in range(len(dn_mutations)):  # iterating through all DN combinations
    substi0 = list(DNA_Nucleotides_set - {wild_type[dn_mutations[dn_comb][0]]})  # all substitutions
    substi1 = list(DNA_Nucleotides_set - {wild_type[dn_mutations[dn_comb][1]]})  # all substitutions

    copy_item = list(wild_type[:])
    for s0 in range(len(substi0)):  # iterating through all substitutions
        for s1 in range(len(substi1)):  # iterating through all substitutions
            counterMutationsDN += 1
            copy_item[dn_mutations[dn_comb][0]] = substi0[s0]
            copy_item[dn_mutations[dn_comb][1]] = substi1[s1]
            NewAA = AASequence(copy_item)  # mutant AA

            diff_count = differing_characters(AA, NewAA)
            if diff_count > 1:
                counterMutationsDN_twoAAmutations += 1
                NewTwoAASetDN.add(NewAA)
            if diff_count == 1:
                NewOneAASetDN.add(NewAA)
            NewAASetDN.add(NewAA)
            if AA == NewAA:   # Check for synonymous mutations
                counterDNsyn += 1
counterDN = len(NewAASetDN) - 1  # Count unique amino acid sequences due to DN mutations
counterDNMore = len(NewAASetDN.difference(NewAASetSN))  # Amino acids accessible only by DN and not SN
TWOAAcounterDNMore = len(NewTwoAASetDN.difference(NewAASetSN))  # Amino acids accessible only by DN and not SN
OneAAcounterDNMore = len(NewOneAASetDN.difference(NewAASetSN))  # Amino acids accessible only by DN and not SN



print("Accessible Amino acids through DN mutations", counterDN)
print("Amino acids accessible only through DN mutations", counterDNMore)

print("One Amino acids accessible only through DN mutations", OneAAcounterDNMore)
print("Two Amino acids accessible only through DN mutations", TWOAAcounterDNMore)



print("\n")
print("Total distinct SN mutations", counterMutationsSN)
print("Expectation", Codons*3*3)
print("\n")
print("Total distinct DN mutations", counterMutationsDN)
print("Expectation", (Codons*3-1)*3*3)

