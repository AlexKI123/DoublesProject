import itertools
from itertools import product
import timeit

# initial and final stop codons are ignored

# Lists of DNA nucleotides and their respective sets.
DNA_Nucleotides = ['A', 'C', 'G', 'T']
DNA_Nucleotides_set = {'A', 'C', 'G', 'T'}

# Dictionary for DNA nucleotide complement mapping.
DNA_ReverseComplement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

# Dictionary of DNA codons to corresponding amino acids.
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


# Function to convert DNA sequence (ntSeq) into amino acid sequence.
def AASequence(ntSeq):
    Codons = int(len(ntSeq)/3)
    AA = []
    for c in range(Codons):
        codonSeq = ''.join(ntSeq[0+3*c:3+3*c])
        AA.append(DNA_Codons[codonSeq])
    return ''.join(AA)

# Fixed number of codons to consider.
Codons = 3
print("RingCodons", Codons)

# Generate all possible DN combinations for the given number of codons.
dn_mutations = [[i,i+1] for i in range(Codons*3-1)]
dn_mutations.append([0, Codons*3-1])

# Print DN combinations and total number of possible genotypes.
print(dn_mutations)
print("Genotypes", 4**(Codons*3))

# Initializing counters to track various mutation outcomes.
counterSN = 0  # counter for AA accessibility through SN
counterDN = 0  # counter for AA accessibility through DN
counterDNMore = 0  # counter for AA accessibility through DN / SN
counterSNMore = 0  # counter for AA accessibility through SN / DN
counterSNandDN = 0  # counter for AA accessibility through SN and DN
counterSNsyn = 0  # counter for synonymous SN mutation
counterDNsyn = 0  # counter for synonymous DN mutation
counter_all = 0
counter_all_SN = 0
counter_all_DN = 0


# Iterating through all possible DNA sequences of the given codon length.
start = timeit.default_timer()
ci = 0
for item in product(DNA_Nucleotides, repeat=3*Codons):
    ci += 1
    if ci % 1_000_000 == 0:
        print(ci)
        stop = timeit.default_timer()
        print(stop-start)

    AA = AASequence(item)  # wild type AA
    if AA.find('_') == -1:  # Ignore sequences with stop codons.
        counter_all += 1
        NewAASetSN = {AA}
        for locus in range(len(item)):  # iterating through all loci
            substi = list(DNA_Nucleotides_set-{item[locus]})  # all substitutions
            copy_item = list(item[:])
            for s in range(len(substi)):  # iterating through all substitutions
                counter_all_SN += 1
                copy_item[locus] = substi[s]
                NewAA = AASequence(copy_item)  # mutant AA
                if NewAA.find('_') == -1:
                    NewAASetSN.add(NewAA)
                    if AA == NewAA:
                        counterSNsyn += 1
        counterSN += (len(NewAASetSN)-1)/Codons  # how many distinct mutant AA are accessible?

        NewAASetDN = {AA}  # wild type AA
        for dn_comb in range(len(dn_mutations)):  # iterating through all DN combinations
            substi0 = list(DNA_Nucleotides_set-{item[dn_mutations[dn_comb][0]]})  # all substitutions
            substi1 = list(DNA_Nucleotides_set-{item[dn_mutations[dn_comb][1]]})  # all substitutions

            copy_item = list(item).copy()
            for s0 in range(len(substi0)):  # iterating through all substitutions
                for s1 in range(len(substi1)):  # iterating through all substitutions
                    counter_all_DN += 1
                    copy_item[dn_mutations[dn_comb][0]] = substi0[s0]
                    copy_item[dn_mutations[dn_comb][1]] = substi1[s1]
                    NewAA = AASequence(copy_item)  # mutant AA
                    if NewAA.find('_') == -1:
                        NewAASetDN.add(NewAA)
                        if AA == NewAA:
                            counterDNsyn += 1
        counterDN += (len(NewAASetDN)-1)/Codons
        counterDNMore += (len(NewAASetDN.difference(NewAASetSN)))/Codons  # how many distinct mutant AA are accessible?
        counterSNMore += (len(NewAASetSN.difference(NewAASetDN)))/Codons  # how many distinct mutant AA are accessible?

        counterSNandDN += (len(NewAASetSN.intersection(NewAASetDN)))/Codons  # how many distinct mutant AA are accessible?

print("4**(Codons*3)", 4**(Codons*3))
print("counter", counter_all)
print("\n")
print("4**(Codons*3)*Codons*3*3", 4**(Codons*3)*Codons*3*3)
print("counter_all_SN", counter_all_SN)
print("\n")
print("4**(Codons*3)*(Codons*3)*3*3", 4**(Codons*3)*(Codons*3)*3*3)
print("counter_all_DN", counter_all_DN)
print("\n")

print("1: ", counterSN/counter_all)  # Acc AA per nucleotide triplet through SN
print("2: ", counterDN/counter_all)  # Acc AA per nucleotide triplet through DN
print("3: ", counterDNMore/counter_all)  # Acc AA per nucleotide triplet through DN / SN
print("4: ", counterSNMore/counter_all)  # Acc AA per nucleotide triplet through SN / DN
print("5: ", counterSNandDN/counter_all)  # Acc AA per nucleotide triplet through DN and SN

print("\n")
print(counterSNsyn/counter_all_SN)  # Prob SN is synonymous
print(counterDNsyn/counter_all_DN)  # Prob DN is synonymous


