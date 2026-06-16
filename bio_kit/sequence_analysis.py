import re
from .data.codon_table import STANDARD_CODON_TABLE, START_CODONS, STOP_CODONS, COMPLEMENT_TABLE
from .sequence_io import Sequence

def complement(seq):
    if isinstance(seq, Sequence):
        seq = seq.seq
    return seq.translate(COMPLEMENT_TABLE)

def reverse_complement(seq):
    if isinstance(seq, Sequence):
        seq = seq.seq
    return complement(seq)[::-1]

def transcribe(dna):
    if isinstance(dna, Sequence):
        dna = dna.seq
    rna = dna.replace('T', 'U')
    return rna

def _translate_codon(codon):
    codon = codon.upper().replace('U', 'T')
    if len(codon) < 3:
        return ''
    if 'N' in codon:
        return 'X'
    return STANDARD_CODON_TABLE.get(codon, 'X')

def translate(seq, frame=0, to_stop=False):
    if isinstance(seq, Sequence):
        seq = seq.seq
    seq = seq.upper().replace('U', 'T')
    if frame < 0:
        seq = reverse_complement(seq)
        frame = abs(frame) - 1
    seq = seq[frame:]
    protein = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if len(codon) < 3:
            break
        aa = _translate_codon(codon)
        if to_stop and aa == '*':
            break
        protein.append(aa)
    return ''.join(protein)

def six_frame_translation(seq):
    if isinstance(seq, Sequence):
        seq = seq.seq
    frames = {}
    for i in range(3):
        frames[f'+{i+1}'] = translate(seq, frame=i)
    rc = reverse_complement(seq)
    for i in range(3):
        frames[f'-{i+1}'] = translate(rc, frame=i)
    return frames

def _find_orfs_frame(seq, frame, min_length=30, start_codons=None, stop_codons=None):
    if start_codons is None:
        start_codons = START_CODONS
    if stop_codons is None:
        stop_codons = STOP_CODONS
    seq = seq.upper().replace('U', 'T')
    orfs = []
    i = frame
    start = -1
    while i < len(seq) - 2:
        codon = seq[i:i+3]
        if start == -1 and codon in start_codons:
            start = i
        elif start != -1 and codon in stop_codons:
            end = i + 3
            length_aa = (end - start) // 3
            if length_aa * 3 >= min_length * 3:
                orfs.append({
                    'start': start,
                    'end': end,
                    'frame': frame,
                    'strand': '+' if frame >= 0 else '-',
                    'length_bp': end - start,
                    'length_aa': length_aa,
                    'protein': translate(seq[start:end])
                })
            start = -1
        i += 3
    return orfs

def find_orfs(seq, min_length=30, start_codons=None, stop_codons=None):
    if isinstance(seq, Sequence):
        seq = seq.seq
    all_orfs = []
    for frame in range(3):
        all_orfs.extend(_find_orfs_frame(seq, frame, min_length, start_codons, stop_codons))
    rc = reverse_complement(seq)
    seq_len = len(seq)
    for frame in range(3):
        orfs = _find_orfs_frame(rc, frame, min_length, start_codons, stop_codons)
        for orf in orfs:
            orf['strand'] = '-'
            orf['frame'] = -(frame + 1)
            start = seq_len - orf['end']
            end = seq_len - orf['start']
            orf['start'] = start
            orf['end'] = end
            all_orfs.append(orf)
    all_orfs.sort(key=lambda x: x['length_aa'], reverse=True)
    return all_orfs
