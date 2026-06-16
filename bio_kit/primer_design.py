import math
from .sequence_analysis import reverse_complement

def calculate_gc(seq):
    if not seq:
        return 0.0
    seq = seq.upper()
    gc = seq.count('G') + seq.count('C')
    return gc / len(seq) * 100

def calculate_tm(primer):
    if not primer:
        return 0.0
    seq = primer.upper()
    if len(seq) <= 14:
        return (seq.count('A') + seq.count('T')) * 2 + (seq.count('G') + seq.count('C')) * 4
    else:
        gc = calculate_gc(seq)
        return 64.9 + 41 * (seq.count('G') + seq.count('C') - 16.4) / len(seq)

def check_dimer(primer, min_bp=4):
    seq = primer.upper()
    rc = reverse_complement(primer)
    max_match = 0
    for i in range(len(seq)):
        for j in range(len(rc)):
            match_len = 0
            k = 0
            while i + k < len(seq) and j + k < len(rc) and seq[i+k] == rc[j+k]:
                match_len += 1
                k += 1
            if match_len > max_match:
                max_match = match_len
    return {'has_dimer': max_match >= min_bp, 'max_match': max_match}

def check_hairpin(primer, min_loop=3, min_stem=4):
    seq = primer.upper()
    n = len(seq)
    max_stem = 0
    best_loop = 0
    
    for stem_len in range(min_stem, n // 2):
        for i in range(n - 2 * stem_len - min_loop + 1):
            loop_start = i + stem_len
            loop_end = n - stem_len
            if loop_end <= loop_start:
                continue
            loop_size = loop_end - loop_start
            if loop_size < min_loop:
                continue
            left = seq[i:i+stem_len]
            right = reverse_complement(seq[n-stem_len-i:n-i])
            if left == right:
                if stem_len > max_stem:
                    max_stem = stem_len
                    best_loop = loop_size
    
    return {
        'has_hairpin': max_stem >= min_stem,
        'max_stem': max_stem,
        'loop_size': best_loop,
    }

def _generate_primers(template, start, end, length_range=(18, 25)):
    forward_primers = []
    min_len, max_len = length_range
    region = template[start:end]
    
    for length in range(min_len, max_len + 1):
        for i in range(len(region) - length + 1):
            primer_seq = region[i:i+length]
            forward_primers.append({
                'sequence': primer_seq,
                'position': start + i,
                'direction': 'forward',
            })
    
    reverse_primers = []
    rc_region = reverse_complement(region)
    for length in range(min_len, max_len + 1):
        for i in range(len(rc_region) - length + 1):
            primer_seq = rc_region[i:i+length]
            orig_pos = end - i - length
            reverse_primers.append({
                'sequence': primer_seq,
                'position': orig_pos,
                'direction': 'reverse',
            })
    
    return forward_primers, reverse_primers

def design_primers(template, start=None, end=None, length_range=(18, 25),
                   tm_range=(55, 65), gc_range=(40, 60),
                   product_size_range=(100, 500),
                   min_product_size=None, max_product_size=None,
                   max_dimers=4, max_hairpin=3):
    if start is None:
        start = 0
    if end is None:
        end = len(template)
    
    if min_product_size is not None:
        product_size_range = (min_product_size, product_size_range[1])
    if max_product_size is not None:
        product_size_range = (product_size_range[0], max_product_size)
    
    template = template.upper()
    
    forward_primers, reverse_primers = _generate_primers(template, start, end, length_range)
    
    scored_forward = []
    for primer in forward_primers:
        seq = primer['sequence']
        tm = calculate_tm(seq)
        gc = calculate_gc(seq)
        dimer = check_dimer(seq, max_dimers)
        hairpin = check_hairpin(seq, min_stem=max_hairpin)
        
        if (tm_range[0] <= tm <= tm_range[1] and 
            gc_range[0] <= gc <= gc_range[1] and
            not dimer['has_dimer'] and
            not hairpin['has_hairpin']):
            
            tm_optimal = (tm_range[0] + tm_range[1]) / 2
            gc_optimal = (gc_range[0] + gc_range[1]) / 2
            score = abs(tm - tm_optimal) + abs(gc - gc_optimal) * 0.5
            
            primer.update({
                'tm': tm,
                'gc': gc,
                'dimer_max': dimer['max_match'],
                'hairpin_max': hairpin['max_stem'],
                'score': score,
            })
            scored_forward.append(primer)
    
    scored_reverse = []
    for primer in reverse_primers:
        seq = primer['sequence']
        tm = calculate_tm(seq)
        gc = calculate_gc(seq)
        dimer = check_dimer(seq, max_dimers)
        hairpin = check_hairpin(seq, min_stem=max_hairpin)
        
        if (tm_range[0] <= tm <= tm_range[1] and 
            gc_range[0] <= gc <= gc_range[1] and
            not dimer['has_dimer'] and
            not hairpin['has_hairpin']):
            
            tm_optimal = (tm_range[0] + tm_range[1]) / 2
            gc_optimal = (gc_range[0] + gc_range[1]) / 2
            score = abs(tm - tm_optimal) + abs(gc - gc_optimal) * 0.5
            
            primer.update({
                'tm': tm,
                'gc': gc,
                'dimer_max': dimer['max_match'],
                'hairpin_max': hairpin['max_stem'],
                'score': score,
            })
            scored_reverse.append(primer)
    
    scored_forward.sort(key=lambda x: x['score'])
    scored_reverse.sort(key=lambda x: x['score'])
    
    pairs = []
    for fp in scored_forward[:20]:
        for rp in scored_reverse[:20]:
            product_size = rp['position'] + len(rp['sequence']) - fp['position']
            if product_size_range[0] <= product_size <= product_size_range[1]:
                tm_diff = abs(fp['tm'] - rp['tm'])
                pair_score = fp['score'] + rp['score'] + tm_diff * 2
                
                cross_dimer = _cross_dimer(fp['sequence'], rp['sequence'])
                if cross_dimer < 5:
                    pairs.append({
                        'forward': fp,
                        'reverse': rp,
                        'product_size': product_size,
                        'tm_diff': tm_diff,
                        'cross_dimer': cross_dimer,
                        'score': pair_score,
                    })
    
    pairs.sort(key=lambda x: x['score'])
    
    return {
        'forward_primers': scored_forward[:10],
        'reverse_primers': scored_reverse[:10],
        'best_pairs': pairs[:10],
    }

def _cross_dimer(seq1, seq2):
    rc2 = reverse_complement(seq2)
    max_match = 0
    for i in range(len(seq1)):
        for j in range(len(rc2)):
            match_len = 0
            k = 0
            while i + k < len(seq1) and j + k < len(rc2) and seq1[i+k] == rc2[j+k]:
                match_len += 1
                k += 1
            if match_len > max_match:
                max_match = match_len
    return max_match
