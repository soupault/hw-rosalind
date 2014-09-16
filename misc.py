# Various utilities
from pprint import pprint
from random import choice


def read_fasta(filename):
    """
    Read file in FASTA format.
    :param filename: full path to FASTA file
    :return: dictionary with labels as 'keys' and genetic strings as 'values'
    """
    with open(filename) as f:
        lines = [e.strip() for e in f.readlines()]

    dnas = dict()
    label_idxs = [idx for idx, line in enumerate(lines) if line.startswith('>')]
    for n in range(len(label_idxs)):
        slice_l = label_idxs[n] + 1
        if n == len(label_idxs) - 1:
            slice_r = len(lines)
        else:
            slice_r = label_idxs[n + 1]

        dnas[lines[label_idxs[n]][1:]] = ''.join(lines[slice_l:slice_r])
    return dnas


def monoisotopic_mass_table():
    mass_table = {'A': 71.03711 , 'C': 103.00919, 'D': 115.02694,
                  'E': 129.04259, 'F': 147.06841, 'G': 57.02146,
                  'H': 137.05891, 'I': 113.08406, 'K': 128.09496,
                  'L': 113.08406, 'M': 131.04049, 'N': 114.04293,
                  'P': 97.05276,  'Q': 128.05858, 'R': 156.10111,
                  'S': 87.03203,  'T': 101.04768, 'V': 99.06841,
                  'W': 186.07931, 'Y': 163.06333}
    return mass_table


def rna_codon_table():
    """
    Function contains matching tables between codons and amino acids
    :return: encoding dictionary {'codon': 'acid/stop', ...}
    """
    code_table = {'GUC': 'V', 'ACC': 'T', 'GUA': 'V', 'GUG': 'V', 'GUU': 'V',
                  'AAC': 'N', 'CCU': 'P', 'UGG': 'W', 'AGC': 'S', 'AUC': 'I',
                  'CAU': 'H', 'AAU': 'N', 'AGU': 'S', 'ACU': 'T', 'CAC': 'H',
                  'ACG': 'T', 'CCG': 'P', 'CCA': 'P', 'ACA': 'T', 'CCC': 'P',
                  'UGU': 'C', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 'UGC': 'C',
                  'CAG': 'Q', 'GAU': 'D', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S',
                  'AGG': 'R', 'GGG': 'G', 'UCC': 'S', 'UCA': 'S', 'UAA': 'Stop',
                  'GGA': 'G', 'UAC': 'Y', 'GAC': 'D', 'GAA': 'E', 'AUA': 'I',
                  'GCA': 'A', 'CUU': 'L', 'GGC': 'G', 'AUG': 'M', 'UGA': 'Stop',
                  'CUG': 'L', 'GAG': 'E', 'CUC': 'L', 'AGA': 'R', 'CUA': 'L',
                  'GCC': 'A', 'AAA': 'K', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F',
                  'CGU': 'R', 'CGA': 'R', 'GCU': 'A', 'UAG': 'Stop', 'AUU': 'I',
                  'UUG': 'L', 'UUA': 'L', 'CGC': 'R', 'UUC': 'F'}

    return code_table


def dna_codon_table():
    code_table = {'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAT': 'N', 'ACA': 'T', 'ACC': 'T',
                  'ACG': 'T', 'ACT': 'T', 'AGA': 'R', 'AGC': 'S', 'AGG': 'R', 'AGT': 'S',
                  'ATA': 'I', 'ATC': 'I', 'ATG': 'M', 'ATT': 'I', 'CAA': 'Q', 'CAC': 'H',
                  'CAG': 'Q', 'CAT': 'H', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
                  'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R', 'CTA': 'L', 'CTC': 'L',
                  'CTG': 'L', 'CTT': 'L', 'GAA': 'E', 'GAC': 'D', 'GAG': 'E', 'GAT': 'D',
                  'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A', 'GGA': 'G', 'GGC': 'G',
                  'GGG': 'G', 'GGT': 'G', 'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
                  'TAA': 'Stop', 'TAC': 'Y', 'TAG': 'Stop', 'TAT': 'Y', 'TCA': 'S', 'TCC': 'S',
                  'TCG': 'S', 'TCT': 'S', 'TGA': 'Stop', 'TGC': 'C', 'TGG': 'W', 'TGT': 'C',
                  'TTA': 'L', 'TTC': 'F', 'TTG': 'L', 'TTT': 'F'}
    return code_table


def reverse(s):
    return s[::-1]


def complement(s, acid='rna'):
    result = ''
    if acid == 'rna':
        table = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}
    elif acid == 'dna':
        table = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    for e in s:
        result += table[e]
    return result


def reverse_complement(s, acid='rna'):
    return complement(reverse(s), acid=acid)


def is_transition(c1, c2):
    if {c1, c2} == {'A', 'G'} or {c1, c2} == {'C', 'T'}:
        return True
    else:
        return False


def is_transversion(c1, c2):
    if {c1, c2} == {'A', 'C'} or {c1, c2} == {'G', 'C'} or \
       {c1, c2} == {'A', 'T'} or {c1, c2} == {'G', 'T'}:
        return True
    else:
        return False


def find_orf(s, acid='rna'):
    if acid == 'rna':
        table = rna_codon_table()
        codon_start = 'AUG'
    elif acid == 'dna':
        table = dna_codon_table()
        codon_start = 'ATG'
    else:
        raise ValueError

    # Find all 'Start' codons
    starts = []
    for idx in range(0, len(s), 3):
        if s[idx:idx + 3] == codon_start:
            starts.append(idx)
    if not starts:
        print("Can't find start codons!")
        return []

    substrings = []
    # For every 'Start' codon
    for idx in range(len(starts)):
        bound = starts[-1] + ((len(s) - starts[idx]) // 3) * 3
        interval = [starts[idx], bound]
        # Look for 'Stop' codon
        for offset in range(interval[0], interval[1], 3):
            if table[s[offset:offset + 3]] == 'Stop':
                substrings.append(s[interval[0]:offset + 3])
                break
    return substrings


def encode_in_protein(s, acid='rna'):
    if acid == 'rna':
        table = rna_codon_table()
    elif acid == 'dna':
        table = dna_codon_table()
    else:
        raise ValueError

    result = ''
    for i in range(0, len(s), 3):
        codon = s[i:i + 3]
        if table[codon] == 'Stop':
            break
        else:
            result += table[codon]
    return result


def divide_array(array, v):
    array_l = [e for e in array if e < v]
    array_v = [e for e in array if e == v]
    array_r = [e for e in array if e > v]
    return array_l, array_v, array_r


def selection(array, k):
    v = choice(array)
    a_l, a_v, a_r = divide_array(array, v)
    if k <= len(a_l):
        return selection(a_l, k)
    elif k > len(a_l) + len(a_v):
        return selection(a_r, k - len(a_l) - len(a_v))
    else:
        return v


def depth_first_search(edges, verts):
    def explore(v):
        visited[verts.index(v)] = True
        # PREVISIT(V)
        for e in edges:
            if v in e:
                u = e[1 - e.index(v)]
                if not visited[verts.index(u)]:
                    explore(u)
        # POSTVISIT(V)

    visited = [False] * len(verts)
    count = 0
    for vert in verts:
        if not visited[verts.index(vert)]:
            explore(vert)
            count += 1
    print("Groups:", count)


def check_bipartiteness(verts, edges, first_elem=1):
    visited = [None] * len(verts)
    colors = [None] * len(verts)
    q = [first_elem]
    colors[verts.index(first_elem)] = True
    while q:
        u = q.pop(0)
        for edge in edges:
            if u in edge:
                v = edge[1 - edge.index(u)]
                if colors[verts.index(v)] is not None and \
                   colors[verts.index(u)] is not None and \
                   colors[verts.index(v)] == colors[verts.index(u)]:
                    return False
                if visited[verts.index(v)] is None:
                    q.append(v)
                    visited[verts.index(v)] = True
                    colors[verts.index(v)] = not colors[verts.index(u)]
    return True


# if __name__ == "__main__":
#     t = None
#     t = t.strip().split()
#     print(t)
#     bob = {}
#     for k, v in zip(t[0::2], t[1::2]):
#         bob[k] = v
#     pprint(bob)
