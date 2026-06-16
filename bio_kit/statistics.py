import math
from collections import Counter
from .data.codon_table import STANDARD_CODON_TABLE

def gc_content(seq):
    if hasattr(seq, 'seq'):
        seq = seq.seq
    seq = seq.upper()
    total = len(seq)
    if total == 0:
        return 0.0
    gc = seq.count('G') + seq.count('C')
    return gc / total * 100

def gc_content_sliding(seq, window=100, step=1):
    if hasattr(seq, 'seq'):
        seq = seq.seq
    seq = seq.upper()
    n = len(seq)
    values = []
    positions = []
    for i in range(0, n - window + 1, step):
        window_seq = seq[i:i+window]
        gc = (window_seq.count('G') + window_seq.count('C')) / window * 100
        values.append(gc)
        positions.append(i + window // 2)
    return positions, values

def n50(sequences):
    lengths = []
    for seq in sequences:
        if hasattr(seq, 'seq'):
            lengths.append(len(seq.seq))
        else:
            lengths.append(len(seq))
    
    lengths.sort(reverse=True)
    total = sum(lengths)
    half = total / 2
    
    cumulative = 0
    for length in lengths:
        cumulative += length
        if cumulative >= half:
            return length
    return 0

def kmer_frequency(seq, k=3):
    if hasattr(seq, 'seq'):
        seq = seq.seq
    seq = seq.upper()
    n = len(seq)
    counts = Counter()
    for i in range(n - k + 1):
        kmer = seq[i:i+k]
        if 'N' not in kmer:
            counts[kmer] += 1
    return dict(counts)

def shannon_entropy(seq):
    if hasattr(seq, 'seq'):
        seq = seq.seq
    seq = seq.upper()
    n = len(seq)
    if n == 0:
        return 0.0
    
    counts = Counter(seq)
    entropy = 0.0
    for count in counts.values():
        p = count / n
        if p > 0:
            entropy -= p * math.log2(p)
    return entropy

def rscu_analysis(seq):
    if hasattr(seq, 'seq'):
        seq = seq.seq
    seq = seq.upper().replace('U', 'T')
    
    codon_counts = Counter()
    amino_counts = Counter()
    
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if len(codon) == 3 and 'N' not in codon:
            aa = STANDARD_CODON_TABLE.get(codon, 'X')
            if aa != 'X':
                codon_counts[codon] += 1
                amino_counts[aa] += 1
    
    rscu = {}
    for codon, count in codon_counts.items():
        aa = STANDARD_CODON_TABLE.get(codon, 'X')
        aa_count = amino_counts.get(aa, 0)
        synonymous = sum(1 for c, a in STANDARD_CODON_TABLE.items() if a == aa)
        if aa_count > 0 and synonymous > 0:
            rscu[codon] = (count / aa_count) * synonymous
        else:
            rscu[codon] = 0.0
    
    return {
        'codon_counts': dict(codon_counts),
        'amino_acid_counts': dict(amino_counts),
        'rscu': rscu,
    }

def sequence_stats(seq):
    if hasattr(seq, 'seq'):
        seq_obj = seq
        seq = seq.seq
    else:
        seq_obj = None
    
    seq = seq.upper()
    
    stats = {
        'length': len(seq),
        'gc_content': gc_content(seq),
        'a_count': seq.count('A'),
        't_count': seq.count('T'),
        'g_count': seq.count('G'),
        'c_count': seq.count('C'),
        'n_count': seq.count('N'),
        'other_count': len(seq) - seq.count('A') - seq.count('T') - seq.count('G') - seq.count('C') - seq.count('N'),
        'shannon_entropy': shannon_entropy(seq),
        'kmer_2': kmer_frequency(seq, k=2),
        'kmer_3': kmer_frequency(seq, k=3),
    }
    
    if seq_obj is not None and seq_obj.quality:
        avg_qual = sum(seq_obj.quality) / len(seq_obj.quality)
        stats['avg_quality'] = avg_qual
    
    return stats
