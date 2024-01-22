import itertools
from itertools import product
import timeit
import random
import math
import numpy as np

# Code to validate the transformations from the nucleotide level to the amino acid level

def string_difference(old_str1, new_str2):
    list1 = list(old_str1)
    list2 = list(new_str2)
    result = [list2[i] for i in range(len(list2)) if list2[i] != list1[i]]
    result=''.join(result)
    result2 = [i for i in range(len(list2)) if list2[i] != list1[i]]
    return result, result2


DNA_Nucleotides = ['A', 'C', 'G', 'T']
DNA_Nucleotides_set = {'A', 'C', 'G', 'T'}

DNA_ReverseComplement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

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

DNA_Subs = {
    "A": ["C", "G", "T"],
    "C": ["A", "G", "T"],
    "G": ["A", "C", "T"],
    "T": ["A", "C", "G"],
}


# computes AA sequences from NT sequence
def AASequence(ntSeq):
    Codons = int(len(ntSeq)/3)
    AA = []
    for c in range(Codons):
        codonSeq = ''.join(ntSeq[0+3*c:3+3*c])
        AA.append(DNA_Codons[codonSeq])
    return ''.join(AA)


'''
########################################################
# # Verify (tilde_nu_s -> tilde_mu_S) transformation # #
########################################################
generations = 10000000
codons = 30
L = codons*3
tilde_nu_S = 10**(-5)
print("U", L*tilde_nu_S)

NTSeq = [random.choice(DNA_Nucleotides) for n in range(L)]
counterTest = 0
counter = 0
for i in range(generations):
    for l in range(L):
        if random.random() < tilde_nu_S:
            counterTest +=1
            previousAA = AASequence(NTSeq)
            NTSeq[l] = random.choice(DNA_Subs[NTSeq[l]])
            laterAA = AASequence(NTSeq)
            if previousAA != laterAA:
                counter +=1


print("tilde_nu_S", tilde_nu_S)
print("Measured tilde_nu_S", counterTest/(generations*L))

print("Measured tilde_mu_S", counter/(generations*codons))
print("Theory tilde_mu_S", tilde_nu_S*3*(1-0.23958333333333334))
print(counter)
'''


'''
#################################################################################
# # Verify (tilde_nu_S -> tilde_mu_S) and (tilde_nu_S -> mu_S) transformation # #
#################################################################################
print("Verify nu_S")
generations = 10000000
print("generations", generations)
codons = 30
L = codons*3
tilde_nu_S = 10**(-3)
print("U", L*tilde_nu_S)

aminoAcidSubDict = dict()

initialcodonset = set()


NTSeq = [random.choice(DNA_Nucleotides) for n in range(L)]
counterTest = 0
counter = 0
for i in range(generations):
    for l in range(L):
        if random.random() < tilde_nu_S:
            counterTest +=1
            previousAA = AASequence(NTSeq)
            previousNT = NTSeq[:]
            NTSeq[l] = random.choice(DNA_Subs[NTSeq[l]])
            laterAA = AASequence(NTSeq)
            if previousAA != laterAA:
                counter += 1
                #print("previousAA", previousAA)
                #print("laterAA", laterAA)
                diffkeyNew, changePosNew = string_difference(previousAA, laterAA)
                diffkeyOld, changePosOld = string_difference(laterAA, previousAA)
                codonnumber = changePosNew[0]
                #print("codonnumber", codonnumber)
                codonblub = previousNT[codonnumber*3] + previousNT[codonnumber*3+1] + previousNT[codonnumber*3+2]
                initialcodonset.add(codonblub)
                #print("codonblub", codonblub)
                #print("previousNT", previousNT)
                #print("previousN2", NTSeq)
                #print("\n")
                diffkey = codonblub + "->" + diffkeyNew
                #print("diffkey", diffkey)
                if diffkey in aminoAcidSubDict.keys():
                    aminoAcidSubDict[diffkey] += 1
                else:
                    aminoAcidSubDict[diffkey] = 1

print("len(initialcodonset)", len(initialcodonset))


print("tilde_nu_S", tilde_nu_S)
print("Measured tilde_nu_S", counterTest/(generations*L))

print("Measured tilde_mu_S", counter/(generations*codons))
print("Theory tilde_mu_S", tilde_nu_S*3*(1-0.23958333333333334))


print(aminoAcidSubDict)
print("Measured mu_S", 4*4*4*np.mean(np.array(list(aminoAcidSubDict.values())))/(generations*codons))
print("Theory mu_S", tilde_nu_S*3*(1-0.23958333333333334)/6.11)

print("ratio", (tilde_nu_S*3*(1-0.23958333333333334))/(4*4*4*np.mean(np.array(list(aminoAcidSubDict.values())))/(generations*codons)))


#print(counter)
'''


'''
#########################################
# # Verify (tilde_nu_D -> tilde_mu_D) # #
#########################################
generations = 1000000
print("generations", generations)
codons = 30
print("codons", codons)
L = codons*3
tilde_nu_D = 1*10**(-3)
print("mud", tilde_nu_D)
print("U", (L-1)*tilde_nu_D)

NTSeq = [random.choice(DNA_Nucleotides) for n in range(L)]
counterTest = 0
counter = 0
for i in range(generations):
    for l in range(L-1):
        if random.random() < tilde_nu_D:
            counterTest +=1
            previousNT = NTSeq.copy()
            previousAA = AASequence(NTSeq)
            NTSeq[l] = random.choice(DNA_Subs[NTSeq[l]])
            NTSeq[l+1] = random.choice(DNA_Subs[NTSeq[l+1]])
            laterAA = AASequence(NTSeq)
            if previousAA != laterAA:
                counter += 1
                break_out_flag = False
                for l in range(L):
                    previousNT_2 = previousNT.copy()
                    for s in range(3):
                        previousNT_2[l] = DNA_Subs[previousNT[l]][s]
                        if AASequence(previousNT_2) == laterAA:
                            counter -= 1
                            break_out_flag = True
                            break
                    if break_out_flag:
                        break


print("tilde_nu_D", tilde_nu_D)
print("Measured tilde_nu_D", counterTest/(generations*(L-1)))

print("Measured tilde_mu_D", counter/(generations*codons))
print("Theory tilde_mu_D", (L-1)/codons*tilde_nu_D*(1-0.011363636363636364)*0.5185)

'''


'''
##################################
# # Verify tilde_nu_D -> mu_D # #
#################################
generations = 10000000
print('generations', generations)
codons = 40
print("codons", codons)
L = codons*3
tilde_nu_D = 10**(-3)
print("U", (L-1)*tilde_nu_D)


aminoAcidSubDict_1 = dict()
aminoAcidSubDict_2 = dict()

initialcodonset_1 = set()
initialcodonset_2 = set()


NTSeq = [random.choice(DNA_Nucleotides) for n in range(L)]
counterTest = 0
counter = 0

counter_1 = 0
counter_2 = 0


for i in range(generations):
    for l in range(L-1):
        if random.random() < tilde_nu_D:
            counterTest +=1
            previousAA = AASequence(NTSeq)
            previousNT = NTSeq[:]
            NTSeq[l] = random.choice(DNA_Subs[NTSeq[l]])
            NTSeq[l+1] = random.choice(DNA_Subs[NTSeq[l+1]])
            laterAA = AASequence(NTSeq)
            if previousAA != laterAA:
                counter += 1
                SNAcc = False
                break_out_flag = False
                for l2 in range(L):
                    previousNT_2 = previousNT.copy()
                    for s in range(3):
                        previousNT_2[l2] = DNA_Subs[previousNT[l2]][s]
                        if AASequence(previousNT_2) == laterAA:
                            counter -= 1
                            SNAcc = True
                            break_out_flag = True
                            break
                    if break_out_flag:
                        break

                if SNAcc == False:
                    diffkeyNew, changePosNew = string_difference(previousAA, laterAA)
                    diffkeyOld, changePosOld = string_difference(laterAA, previousAA)

                    #if len(changePosOld) == 1:
                    if (l+1)%3 != 0:
                        counter_1 +=1
                        codonnumber = changePosNew[0]
                        codonblub = previousNT[codonnumber * 3] + previousNT[codonnumber * 3 + 1] + previousNT[codonnumber * 3 + 2]

                        initialcodonset_1.add(codonblub)
                        diffkey = codonblub + "->" + diffkeyNew
                        if diffkey in aminoAcidSubDict_1.keys():
                            aminoAcidSubDict_1[diffkey] += 1
                        else:
                            aminoAcidSubDict_1[diffkey] = 1

                    #elif len(changePosOld) == 2:
                    else:
                        counter_2 +=1

                        codonblub = previousNT[l-2] + previousNT[l-1] + previousNT[l] + previousNT[l+1] + previousNT[l+2] + previousNT[l+3]
                        diffAA = laterAA[l//3] + laterAA[l//3+1]
                        initialcodonset_2.add(codonblub)
                        diffkey = codonblub + "->" + diffAA
                        if diffkey in aminoAcidSubDict_2.keys():
                            aminoAcidSubDict_2[diffkey] += 1
                        else:
                            aminoAcidSubDict_2[diffkey] = 1



print("counter_1", counter_1)
print("counter_2", counter_2)

print("len(initialcodonset_1)", len(initialcodonset_1))
print("len(initialcodonset_2)", len(initialcodonset_2))

print("len(aminoAcidSubDict_1)", len(aminoAcidSubDict_1))
print("len(aminoAcidSubDict_2)", len(aminoAcidSubDict_2))

print("tilde_nu_D", tilde_nu_D)
print("Measured tilde_nu_D", counterTest/(generations*(L-1)))


avg1=64*np.mean(np.array(list(aminoAcidSubDict_1.values())))/(generations*codons)
if len(aminoAcidSubDict_2) > 0:
    avg2=64*64*np.mean(np.array(list(aminoAcidSubDict_2.values())))/(generations*(codons-1))
else:
    avg2=0
print("avg1", avg1)
print("avg2", avg2)

print("Measured mu_S", 2*codons/(L-1)*avg1+(codons-1)/(L-1)*avg2)
print("Theory mu_D", (L-1)/codons*tilde_nu_D*(1-0.011363636363636364)*0.5185/10.375)

'''



'''
#######################
### Compute Prob[not a acc. by SN]
#######################
generations = 10000000
print("generations", generations)
codons = 30
print("codons", codons)
L = codons*3
mud = 1*10**(-5)
print("mud", mud)
print("U", (L-1)*mud)

NTSeq = [random.choice(DNA_Nucleotides) for n in range(L)]
counterTest = 0
counter = 0
counter_1 = 0
counter_2 = 0
for i in range(generations):
    for l in range(L-1):
        if random.random() < mud:
            counterTest +=1
            previousNT = NTSeq.copy()
            previousAA = AASequence(NTSeq)
            NTSeq[l] = random.choice(DNA_Subs[NTSeq[l]])
            NTSeq[l+1] = random.choice(DNA_Subs[NTSeq[l+1]])
            laterAA = AASequence(NTSeq)
            if previousAA != laterAA:
                counter += 1
                counter_1 += 1
                break_out_flag = False
                for l in range(L):
                    previousNT_2 = previousNT.copy()
                    for s in range(3):
                        previousNT_2[l] = DNA_Subs[previousNT[l]][s]
                        if AASequence(previousNT_2) == laterAA:
                            counter -= 1
                            counter_2 += 1
                            break_out_flag = True
                            break
                    if break_out_flag:
                        break

print(counter_2/counter_1)
print(1-counter_2/counter_1)
'''
