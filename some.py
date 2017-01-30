"""
This module takes forward nucleotide sequence of double strain DNA
and retuns a list of sequences of single strain DNA subsequences
that constitute original double strain DNA when put together.

The idea is to split double strain DNA into single strain 
semi-overlaping pieces with segments of overlap having reasonably
similar melting temperatures. This way the original DNA sequence can
be recreated from these pieces by mixing them together and iteratively
heating them to a common melting temperature and cooling down.

This method of building DNA sequences may reduce costs of building DNA
sequences de novo.
"""
import sys

import math

from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


def melt_temperature_for_sequence_awful(sequence):
    """Return sequence's melting temperature.

    For sequences of 15 to 45 nucleotides only.
    """
    t = 0
    for nucleotide in sequence:
        t += melt_temperature(nucleotide)
    return t


def melt_temperature_for_sequence(sequence):
    """Return sequence's melting temperature.

    For sequences longer than 50 nucleotides.
    """
    gc = 0
    for nucleotide in sequence:
        if nucleotide == "G" or nucleotide == "C":
            gc +=1

    # Tm = 81.5 + 16.6 log M + 41(XG+XC) - 500/L
    m = 0.9
    
    return 81.5 + 16.6 * math.log(m) + 41*gc - float(500)/len(sequence)


def temperature_range(seq):
    """Choose appropriate range for melting temperatures."""
    ts = temperatures_of_subsequences(seq)
    ts_avg = (sum(ts) / float(len(ts)))
    delta = 1.5
    t_min = ts_avg - delta
    t_max = ts_avg + delta
    return t_min, t_max
  

def temperatures_of_subsequences(sequence):
    """Evaluate melting temperatures for subsequences in sequence.

    Return a list of melting temperatures for consecutive step
    nucleotide long subsequences of the sequence.
    """
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
    """Split sequence into isoterm subsequences.
    
    Return a list of indexes delimiting subsequences of melting
    temperatures within range of t_min to t_max.
    """    
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
                
        # print("found!", sequence[i:j], j-i)
        intervals.append(str(sequence[i:j]))
        i = j

    assert(len("".join(intervals)) == len(sequence))
    # print("intervals:", len(intervals))

    return intervals


def seq_fuse(full):
    """
    Return [subsequence1+subsequence2, subsequence2+subsequence3, ... ].
    """
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
        sub_seq = full[a] + full[b]
        read.append (sub_seq)
        a += 1
        b += 1
    return read


def rev_complement(full):
    """
    Reverse complement even subsequences in an input list.
    """
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


def to_single_strain_subsequences(sequence, acceptable_range=6):
    """
    Split double strain sequence into single strain subsequences of
    melting temperatures that are within temperature_range.
    """
    seq = Seq(sequence)

    t_min, t_max = temperature_range(seq)

    intervals = find_intervals(seq, t_min, t_max)

    # tot_seq = melt_point(subseqs)                               # ['TAACATGGCT','CTGTGGCATG','AAATTTACGG','TATACGCACA','AAAAAATCGA','CCAATCGATA','AGTCATCGGT']

    fus_seq = seq_fuse(intervals)
    com_seq = rev_complement(fus_seq)

    subsequences = [str(patch) for patch in com_seq]

    temps = [mt.Tm_NN(item) for item in subsequences]

    print("Minimal temperature:", min(temps))
    print("Maximal temperature:", max(temps))
    print("Temperature range:", max(temps) - min(temps))

    first_strain = [str(subseq) for index, subseq in enumerate(com_seq) if index % 2 == 0]
    second_strain = [str(subseq) for index, subseq in enumerate(com_seq) if index % 2 == 1]

    # todo make assertions
    print(len("".join(first_strain)), len(str(seq)))
    print(len("".join(second_strain)), len(str(seq)))

    # test if sequence can be constructed 

    return subsequences


def main(argv):
    if argv:
        sequence = argv[0]
    else:
        try:        
            with open('input.txt','r') as f:
                sequence = f.readline().rsplit('\n')[0]
        except (OSError, IOError) as e:
            print("DNA sequence file not found. ", e)
            return

    print("Analysing sequence", sequence[:100], "...")
    subsequences = to_single_strain_subsequences(sequence)
    print("Results:")
    for subsequence in subsequences:
        print(subsequence)
    # for index, subsequence in enumerate(subsequences):
    #     print(index, subsequence)


if __name__ == '__main__':
    main(sys.argv[1:])