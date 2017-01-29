import math

from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

# # paste your code pls
# tr="GATATATGCATATACTT"
# t = "ATAT"

# lenT= len(t)


# for i in range(len(tr)-len(t)+1):
#     print("     ", tr[i])
#     if tr[i:i+lenT] == t:
#         print(i)

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

def seq_fuse (full):
    a = 0
    b = 1
    read = []
    sub_seq = []
    leng = len(full)
    if leng % 2 != 0:
        full[leng - 2] = full[leng - 2] + full[leng - 1]
        del full[leng - 1]
        leng = len(full)
    while b <= leng - 1:
        sub_seq = full [a] + full [b]
        read.append (sub_seq)
        a += 1
        b += 1
    return read

def rev_complement (full):
    i = 0
    sub = []
    res = []
    leng = len(full) - 1
    while i <= leng:
        sub = full[i]
        if i % 2 == 0:
            res.append(sub)
        else:
            sub_seq = Seq(full[i])
            sub = sub_seq.reverse_complement()
            res.append(sub)
        i += 1
    return res



seq = Seq("ATATATATATGACTCCCCGTCCATATGTTACTTACCCCGATGAATATCAATATGTCGCTTAGGTGTGGTCACTCTGATATTATTCATATCATAGACACCACCATCCTGACTGATTGGTTTATCATGATGTAATTCAAAGCGTTCCCTACCACCTACTTGGTCTTTCTTCCTTGCAAAAGGTGCTTTTCCCTTTTGAATGTTCGTCTTATTACTGCCTTTAAATTGCTTACTAAGATCGGGATCTTTTGACACTTCTTCCCAAAACTTCTTCCGGAAATCGTCAAAGTTTTTAAATTCTTTATCACGCAACTTATCAGCAATGCGATCTGGAATTGGCGCTCCTGAATCTTTACCTGCATCATCCAGCCATTTATCACCAACTGGTTTACCTTTACCTGTCGCCTTCCCTGGCTTATTCCGTTTACTCTCCTTATCTAATTTATCCTTAGCGTCCTTTTCTTTATTTTCCTTCTGTTTGCGGCGCTCCTGCGCGGCACTTAATGCAGCATCAGCATCTGACTTCTCTTTTGCAGCAGCATCAAATGCAGCCTGCTTATTATTTACATCCGTCTGCGCCCGCTGGGCTTTAAGCCCGGCCATTTGCCACATTCTGTGACCGCCAGCCATTGGGTCATGGGCAAATCGATCAAATTGTTTTATTTCAGCTATTGCATCAGCAAGAGTTTTATTCGCTGCATCAAGTTCGCTTTTACGCGAATTATAAACCTGAACAGCTTTAGCCTGTCGCTCCTGATTTCTGGCAACATCTTCATTTGCCTGATTCAGCTCTGCACGCGCGCGTTCATAATTTCGCTCAGCCGCTTCAACCGGATGCGTAGCATCCCATTCCTGCTGACGGCGATTTTCTTCATCCTGACGTTGTTTTACCTGGTCAGGGCTAAGAACATCACTCACTGAAACATATACGGCATTATGACCGCTGTCCTTCGGGAATCGAATAACTGCATCCCTGGTATTACCACCCTGAGTAAATCCTGCCGGGCGAACATCCTTATCAGTATTATTTGTAACACCTGGGCTTAATGTCTGTACTGCTGGCGTACTGTTATTAACTGAAATATTCAGAACAGGTGCACCTGGAATTGATGCCGTAAAAACACCTGGACGTTCGGTAGGTTTTGCATCAACCACCGGAACACTCATCGGAACACCTGAAACAACCGAAATATTCTGTCGTTCGTCTTTTACATCATCAACAACACGAACATTTACGTTTACTGTTGCCTTATCGAGAGGTAATGAACTGACAGGTGATTCAGTAATATCATCTGCGGGTAATGACGTCACAATCTTTGACATCATATTGGGGTCATCTTTCGCTATTTGTGATGGCAATACACCATATAAAGCCACCCCCCAAAGACCAAATTTAAACGGTCCTTTCAGGGCAGCCATAATATCAGCAATAGCTGCCGATAATGCTCCCGCTGAAATACTGACCGCCAGACCGCCAGCTCCTGGAGTGGAAAGTGCCGGAAAACCAAATGCCACTGGCGCAGCTACTGCTGACAGATTACCGCCTGTTCCCGAGCCACCACCGGAATTACTATTCCCCCCGCCATTACCATGACCGGAACCACCACCCCAGTGAATGCCGCTACCGGAACCACCACCCCACGGGTTATTTTCCGAACTCCACCCGGAACCATCAGAAGCACCACCACCTACACCAAGCCCGGTCGGGCCACCATTAATGTTACCACTTGTGCTATGCGCGCCCGTGTTATGGCCGCGTCCATCGCCACCGCTCATGGATCCACCTCCTGAATTCCATATGAACGTTGTCACGGAAATAATAGGTTCCCTACTGGCACTCCTGTCAGTACTTCCATAGCATTTTCTTCACTTCCTTACCGCTCGGTTTTTGCCGTTCGGCCCCAAATTTGTCCCGACTGGATCCACCTCCTGAATTCGAAAACGTTCATATATATATA")



ts = print_temperatures(seq)

ts_avg = (sum(ts) / float(len(ts)))
delta = 1.5
t_min = ts_avg - delta
t_max = ts_avg + delta

find_intervals(seq,t_min,t_max)

tot_seq = melt_point(new_seq)                               # ['TAACATGGCT','CTGTGGCATG','AAATTTACGG','TATACGCACA','AAAAAATCGA','CCAATCGATA','AGTCATCGGT']
fus_seq = seq_fuse(tot_seq)
com_seq = rev_complement(fus_seq)



# plt.hist(ts)

# plt.show()
2   