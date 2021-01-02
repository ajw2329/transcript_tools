import sys
import gzip
import random



def import_fasta(fasta_path):
    '''
    Turns a fasta file into a dictionary of sequences
    indexed by sequence name (i.e. chromosome name in the
    case of a genome fasta file)

    Parameters
    ----------
    fasta_path : str
        Path to the fasta file

    Returns
    -------
    fasta_dict : dict
        Dictionary containing sequences indexed by their
        names
    '''

    if fasta_path.endswith(".gz"):

        file = gzip.open(fasta_path, 'rt')

    else:

        file = open(fasta_path, 'r')

    fasta_dict = {}

    for line in file:

        if line.startswith(">"):

            seqname = line[1:].strip().split()[0]

            print(seqname)

            fasta_dict.setdefault(seqname, [])

        else:

            fasta_dict[seqname].append(line.strip())

    file.close()

    for seqname, seq in fasta_dict.items():

        fasta_dict[seqname] = ''.join(seq)

    return fasta_dict


def calc_seqname_len(fasta_dict):

    seqname_len_dict = {}

    for seqname, seq in fasta_dict.items():

        seqname_len_dict[seqname] = len(seq)

    return seqname_len_dict



def slice_seq(seqname, start, end, fasta_dict):
    '''
    Retrieves a subsequence from a dictionary of
    full sequences using the sequence name and the
    1-based coordinates of the sequence.

    Parameters
    ----------
    seqname : str
        The name of the sequence (i.e. chromosome name
        in the case of a genome)
    start : int
        1-based coordinate of the beginning of the
        subsequence
    end : int
        1-based coordinate of the end of the subsequence
    fasta_dict : dict
        Dictionary containing full sequences indexed by name
        (i.e. chromosomes in the case of a genome)

    Returns
    -------
    seq : str
        A string matching the subsequence of interest.

    >>> slice_seq('chrA', 5, 10, {'chrA': 'AGGACCGATATGTTGCACGTTGAC'})
    'CCGATA'
    '''

    seq = fasta_dict[seqname][start - 1: end]

    return seq


def gen_random_seq(seq_length):

    bases = ["A", "T", "G", "C"]

    seq = ''.join([random.choice(bases) for i in range(0, seq_length)])

    return seq

def gen_complement(seq):

    seq_list = list(seq)

    comp_list = []

    bp = {"A": "T", "G": "C", "T": "A", "C": "G"}

    comp = ''.join([bp[i] for i in seq_list])

    return comp



def gen_rev_comp(seq):

    seq_list = list(seq)

    seq_list.reverse()

    bp = {"A": "T", "G": "C", "T": "A", "C": "G"}

    rev_comp = ''.join([bp[i] for i in seq_list])

    return rev_comp


def get_non_overlapping_triplets(seq):
    '''
    Takes a sequence and converts to non-overlapping
    triplets, starting at the beginning of the sequence. 
    Triplets should correspond to codons, if the first 3
    nucleotides in the sequence were to correspond to the
    start codon. The resulting triplets are returned as 
    a list of strings, each with length 3. If the 
    provided sequence has length not divisible by 3,
    the leftover nucleotides will be provided as one
    final element in the list of length 1 or 2 as
    necessary.

    Parameters
    ----------

    seq : str
        A string corresponding to a nucleic acid sequence

    Returns 
    -------

    list
        A list of non-overlapping triplets constituting the provided sequence

    >>> get_non_overlapping_triplets("ATGATGATGATGATGAG")                  
    ['ATG', 'ATG', 'ATG', 'ATG', 'ATG', 'AG']  

    >>> get_non_overlapping_triplets("ATGATGATGATGA") 
    ['ATG', 'ATG', 'ATG', 'ATG', 'A'] 

    >>> get_non_overlapping_triplets("ATGATGATGATGATG") 
    ['ATG', 'ATG', 'ATG', 'ATG', 'ATG']
    '''

    triplets = [seq[i:(i+3)] for i in range(0, len(seq), 3)]

    return triplets

