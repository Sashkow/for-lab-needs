import math

from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq

import matplotlib.pyplot as plt
mystring = 'CGTTCCAAAGATGTGGGCATGAGCTTAC'


# # paste your code pls
# tr="GATATATGCATATACTT"
# t = "ATAT"

# lenT= len(t)


# for i in range(len(tr)-len(t)+1):
#     print("     ", tr[i])
#     if tr[i:i+lenT] == t:
#         print(i)

def melt_temperature(nucleotide):
    if nucleotide == "A" or nucleotide == "T":
        return 2
    else:
        return 4

def melt_temperature_for_sequence_awful(seuqence):
    """
    for sequences of 15 to 45 nucleotides only
    """
    t = 0
    for nucleotide in seuqence:
        t += melt_temperature(nucleotide)
    return t

def melt_temperature_for_sequence(sequence):
    """
    for sequences of 15 to 45 nucleotides only
    """
    gc = 0
    for nucleotide in sequence:
        if nucleotide == "G" or nucleotide == "C":
            gc +=1

    # Tm = 81.5 + 16.6 log M + 41(XG+XC) - 500/L
    m = 0.9
    
    return 81.5 + 16.6 * math.log(m) + 41*gc - float(500)/len(sequence)
    




def print_temperatures(sequence):
    i = 0
    step = 50
    temperatures = []
    while i < len(sequence)-step+1:
        subsequence = sequence[i:i+step]
        t = mt.Tm_NN(subsequence)
        temperatures.append(t)

        # print(
        #         subsequence,
        #         t             
        #     )
        i += step
    return temperatures

def find_intervals(sequence, t_min, t_max):
    print(len(sequence), t_min, t_max)
    intervals = []
    i = 0 
    maxstep = 100
    minstep = 10
    while i < len(sequence):
    #     subsequence = sequence[i:i+step]

    #     t = mt.Tm_NN(subsequence)
    #     j = i + step
    #     less = False
    #     while t > t_max or t < t_min:
    #         # print(i,j)
    #         j+=1
    #         subsequence = sequence[i:j]
    #         t = mt.Tm_NN(subsequence)
    #         if j > len(sequence):
    #             print("less at", i)
    #             less = True
    #             break

        # if i + step > len(sequence):
        #     print(sequence[i:i+step], len(sequence[i:i+step]))
        #     intervals.append(str(sequence[i:i+step]))
        #     break

        subsequence = sequence[i:i+maxstep]
        t = mt.Tm_NN(subsequence)
        j = i + len(subsequence)
        while t > t_max or t < t_min:
            # print(i,j)
            j-=1
            subsequence = sequence[i:j]
            t = mt.Tm_NN(subsequence)
            if j < i + minstep:
                print("halepa at", i)
                return               
                
        print("found!", sequence[i:j], j-i)
        intervals.append(str(sequence[i:j]))
        i = j

    assert(len("".join(intervals)) == len(sequence))

    return intervals

    # sprint("intervals:", len(intervals))


seq = Seq("ATATATATATGACTCCCCGTCCATATGTTACTTACCCCGATGAATATCAATATGTCGCTTAGGTGTGGTCACTCTGATATTATTCATATCATAGACACCACCATCCTGACTGATTGGTTTATCATGATGTAATTCAAAGCGTTCCCTACCACCTACTTGGTCTTTCTTCCTTGCAAAAGGTGCTTTTCCCTTTTGAATGTTCGTCTTATTACTGCCTTTAAATTGCTTACTAAGATCGGGATCTTTTGACACTTCTTCCCAAAACTTCTTCCGGAAATCGTCAAAGTTTTTAAATTCTTTATCACGCAACTTATCAGCAATGCGATCTGGAATTGGCGCTCCTGAATCTTTACCTGCATCATCCAGCCATTTATCACCAACTGGTTTACCTTTACCTGTCGCCTTCCCTGGCTTATTCCGTTTACTCTCCTTATCTAATTTATCCTTAGCGTCCTTTTCTTTATTTTCCTTCTGTTTGCGGCGCTCCTGCGCGGCACTTAATGCAGCATCAGCATCTGACTTCTCTTTTGCAGCAGCATCAAATGCAGCCTGCTTATTATTTACATCCGTCTGCGCCCGCTGGGCTTTAAGCCCGGCCATTTGCCACATTCTGTGACCGCCAGCCATTGGGTCATGGGCAAATCGATCAAATTGTTTTATTTCAGCTATTGCATCAGCAAGAGTTTTATTCGCTGCATCAAGTTCGCTTTTACGCGAATTATAAACCTGAACAGCTTTAGCCTGTCGCTCCTGATTTCTGGCAACATCTTCATTTGCCTGATTCAGCTCTGCACGCGCGCGTTCATAATTTCGCTCAGCCGCTTCAACCGGATGCGTAGCATCCCATTCCTGCTGACGGCGATTTTCTTCATCCTGACGTTGTTTTACCTGGTCAGGGCTAAGAACATCACTCACTGAAACATATACGGCATTATGACCGCTGTCCTTCGGGAATCGAATAACTGCATCCCTGGTATTACCACCCTGAGTAAATCCTGCCGGGCGAACATCCTTATCAGTATTATTTGTAACACCTGGGCTTAATGTCTGTACTGCTGGCGTACTGTTATTAACTGAAATATTCAGAACAGGTGCACCTGGAATTGATGCCGTAAAAACACCTGGACGTTCGGTAGGTTTTGCATCAACCACCGGAACACTCATCGGAACACCTGAAACAACCGAAATATTCTGTCGTTCGTCTTTTACATCATCAACAACACGAACATTTACGTTTACTGTTGCCTTATCGAGAGGTAATGAACTGACAGGTGATTCAGTAATATCATCTGCGGGTAATGACGTCACAATCTTTGACATCATATTGGGGTCATCTTTCGCTATTTGTGATGGCAATACACCATATAAAGCCACCCCCCAAAGACCAAATTTAAACGGTCCTTTCAGGGCAGCCATAATATCAGCAATAGCTGCCGATAATGCTCCCGCTGAAATACTGACCGCCAGACCGCCAGCTCCTGGAGTGGAAAGTGCCGGAAAACCAAATGCCACTGGCGCAGCTACTGCTGACAGATTACCGCCTGTTCCCGAGCCACCACCGGAATTACTATTCCCCCCGCCATTACCATGACCGGAACCACCACCCCAGTGAATGCCGCTACCGGAACCACCACCCCACGGGTTATTTTCCGAACTCCACCCGGAACCATCAGAAGCACCACCACCTACACCAAGCCCGGTCGGGCCACCATTAATGTTACCACTTGTGCTATGCGCGCCCGTGTTATGGCCGCGTCCATCGCCACCGCTCATGGATCCACCTCCTGAATTCCATATGAACGTTGTCACGGAAATAATAGGTTCCCTACTGGCACTCCTGTCAGTACTTCCATAGCATTTTCTTCACTTCCTTACCGCTCGGTTTTTGCCGTTCGGCCCCAAATTTGTCCCGACTGGATCCACCTCCTGAATTCGAAAACGTTCATATATATATA")



ts = print_temperatures(seq)

ts_avg = (sum(ts) / float(len(ts)))
delta = 1.5
t_min = ts_avg - delta
t_max = ts_avg + delta

find_intervals(seq,t_min,t_max)





# plt.hist(ts)

# plt.show()
