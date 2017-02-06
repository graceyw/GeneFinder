# -*- coding: utf-8 -*-
"""
Gene Finder Project

@author: Gracey Wilson

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):


    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """

    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'G':
        return 'C'

def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    seq = ""
    for nucleotide in dna:
        comp = get_complement(nucleotide)
        seq = seq + comp
        index = len(seq) - 1
    seq2 = ""
    while index >= 0:
        nucleotide = seq[index]
        index = index - 1
        seq2 = seq2 + nucleotide
    return seq2

def rest_of_ORF(dna):
    """
    Takes a DNA sequence that is assumed to begin with a start
    codon and returns the sequence up to but not including the
    first in frame stop codon.  If there is no in frame stop codon,
    returns the whole string.

    dna: a DNA sequence
    returns: the open reading frame represented as a string

    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    for i in range(3,len(dna)-2,3):
        codon = dna[i:i+3]
        if codon in ['TAG','TAA','TGA']:
            return dna[:i]
    return dna
    #  if dna[3:6] == 'TAG' or dna[3:6] == 'TAA' or dna[3:6] == 'TGA':
    #      return dna[0:3]
    #  elif dna[6:9] == 'TAG' or dna[6:9] == 'TAA' or dna[6:9] == 'TGA':
    #      return dna[0:6]
    #  elif dna[9:12] == 'TAG' or dna[9:12] == 'TAA' or dna[9:12] == 'TGA':
    #      return dna[0:9]
    #  elif dna[12:15] == 'TAG' or dna[12:15] == 'TAA' or dna[12:15] == 'TGA':
    #      return dna[0:12]
    #  else:
    #      return(dna[0:])

def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
         sequence and returns them as a list.  This function should
         only find ORFs that are in the default frame of the sequence
         (i.e. they start on indices that are multiples of 3).
         By non-nested we mean that if an ORF occurs entirely within
         another ORF, it should not be included in the returned list of ORFs.

         dna: a DNA sequence
         returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("ATGTAG")
    ['ATG']
    """
    all_ORFs = []
    i = 0
    while i <= len(dna)-2:
        codon = dna[i:i+3]
        if codon == 'ATG':
            current_ORF = rest_of_ORF(dna[i:])
            all_ORFs.append(current_ORF)
            i = i + len(current_ORF)
        i += 3 # shorthand for i = i + 3
    return all_ORFs
         #don't stop at return statements - keep going through rest of line

def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    all_ORFs = []
    for i in range(0,3):
        all_ORFs_thisframe = find_all_ORFs_oneframe(dna[i:])
        all_ORFs.extend(all_ORFs_thisframe)
    return all_ORFs

def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    both_strands_all_ORFs = find_all_ORFs(dna)
    both_strands_all_ORFs.extend(find_all_ORFs(get_reverse_complement(dna)))
    return both_strands_all_ORFs

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    n = find_all_ORFs_both_strands(dna)
    return max(n,key=len)

def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    length_longest = 0
    for i in range(num_trials):
        shuffled = shuffle_string(dna)
        ORFs_shuffled = find_all_ORFs_both_strands(shuffled)
        longest = max(ORFs_shuffled,key=len)    # gives longest string itself
        current_length_longest = len(longest)
        if current_length_longest > length_longest:
            length_longest = current_length_longest
    return length_longest

def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
              the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    AAseq = ''
    num_codons = len(dna) - len(dna) % 3 -2
    for i in range(0,num_codons,3):
        codon = dna[i:i+3]
        AA = aa_table[codon]
        AAseq = AAseq + AA
    return AAseq

def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    all_AAs = []
    threshold = longest_ORF_noncoding(dna,1500)
    ORFs = find_all_ORFs_both_strands(dna)
    for ORF in ORFs:
        if len(ORF) > threshold:
            all_AAs.append(coding_strand_to_AA(ORF))
    return all_AAs

if __name__ == "__main__":
    dna = load_seq("./data/X73525.fa")
    answer = gene_finder(dna)
    for i in answer:
        print(i)
        print('_____')
    # import doctest
    # doctest.testmod(verbose=True)
