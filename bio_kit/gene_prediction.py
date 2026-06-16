import re
from .sequence_analysis import reverse_complement

class GenePrediction:
    def __init__(self):
        pass

def _weight_matrix_score(seq, matrix):
    score = 0.0
    for i, base in enumerate(seq):
        if i < len(matrix):
            score += matrix[i].get(base, 0)
    return score

def predict_promoters(seq, threshold=-5.0):
    if hasattr(seq, 'seq'):
        seq = seq.seq
    seq = seq.upper()
    
    tata_matrix = [
        {'T': 1.0, 'A': 0.2, 'G': -1.0, 'C': -1.0},
        {'A': 1.0, 'T': 0.2, 'G': -1.0, 'C': -1.0},
        {'T': 0.8, 'A': 0.5, 'G': -0.5, 'C': -0.5},
        {'A': 1.0, 'T': 0.8, 'G': -1.0, 'C': -1.0},
        {'A': 1.0, 'T': 0.5, 'G': -0.5, 'C': -0.5},
        {'T': 0.8, 'A': 0.3, 'G': -1.0, 'C': -1.0},
    ]
    
    pribnow_matrix = [
        {'T': 1.0, 'A': 0.3, 'G': -0.5, 'C': -0.5},
        {'A': 0.8, 'T': 0.5, 'G': -0.5, 'C': -0.5},
        {'T': 1.0, 'A': 0.5, 'G': -0.5, 'C': -0.5},
        {'A': 1.0, 'T': 0.8, 'G': -0.5, 'C': -0.5},
        {'T': 0.8, 'A': 0.5, 'G': -0.5, 'C': -0.5},
        {'T': 0.5, 'A': 0.5, 'G': -0.3, 'C': -0.3},
    ]
    
    minus35_matrix = [
        {'T': 0.8, 'A': 0.5, 'G': -0.3, 'C': -0.3},
        {'T': 1.0, 'A': 0.3, 'G': -0.5, 'C': -0.5},
        {'G': 0.8, 'T': 0.3, 'A': -0.3, 'C': -0.3},
        {'A': 0.8, 'C': 0.5, 'G': 0.3, 'T': -0.3},
        {'C': 0.8, 'A': 0.5, 'T': -0.3, 'G': -0.3},
        {'A': 0.5, 'T': 0.3, 'G': -0.2, 'C': -0.2},
    ]
    
    promoters = []
    n = len(seq)
    
    for i in range(n - 6):
        window = seq[i:i+6]
        if 'N' in window:
            continue
        
        tata_score = _weight_matrix_score(window, tata_matrix)
        pribnow_score = _weight_matrix_score(window, pribnow_matrix)
        minus35_score = _weight_matrix_score(window, minus35_matrix)
        
        max_score = max(tata_score, pribnow_score, minus35_score)
        
        if max_score >= threshold:
            promoter_type = 'TATA'
            if max_score == pribnow_score:
                promoter_type = 'Pribnow (-10)'
            elif max_score == minus35_score:
                promoter_type = '-35'
            
            promoters.append({
                'start': i,
                'end': i + 6,
                'score': max_score,
                'tata_score': tata_score,
                'pribnow_score': pribnow_score,
                'minus35_score': minus35_score,
                'sequence': window,
                'strand': '+',
                'type': promoter_type,
            })
    
    rc = reverse_complement(seq)
    for i in range(n - 6):
        window = rc[i:i+6]
        if 'N' in window:
            continue
        
        tata_score = _weight_matrix_score(window, tata_matrix)
        pribnow_score = _weight_matrix_score(window, pribnow_matrix)
        minus35_score = _weight_matrix_score(window, minus35_matrix)
        
        max_score = max(tata_score, pribnow_score, minus35_score)
        
        if max_score >= threshold:
            orig_start = n - i - 6
            promoter_type = 'TATA'
            if max_score == pribnow_score:
                promoter_type = 'Pribnow (-10)'
            elif max_score == minus35_score:
                promoter_type = '-35'
            
            promoters.append({
                'start': orig_start,
                'end': orig_start + 6,
                'score': max_score,
                'tata_score': tata_score,
                'pribnow_score': pribnow_score,
                'minus35_score': minus35_score,
                'sequence': window,
                'strand': '-',
                'type': promoter_type,
            })
    
    promoters.sort(key=lambda x: x['score'], reverse=True)
    return promoters

def predict_terminators(seq, stem_length=6, loop_min=3, loop_max=8, u_stretch=4):
    if hasattr(seq, 'seq'):
        seq = seq.seq
    seq = seq.upper()
    
    terminators = []
    n = len(seq)
    
    for i in range(n - 2 * stem_length - loop_min):
        best_terminator = None
        best_stem = 0
        
        for loop_size in range(loop_min, loop_max + 1):
            left_start = i
            left_end = i + stem_length
            right_end = i + stem_length + loop_size + stem_length
            
            if right_end > n:
                break
            
            right_start = i + stem_length + loop_size
            left = seq[left_start:left_end]
            right = seq[right_start:right_end]
            
            rc_right = reverse_complement(right)
            matches = sum(1 for a, b in zip(left, rc_right) if a == b)
            
            if matches >= stem_length * 0.7 and matches > best_stem:
                u_start = right_end
                u_end = u_start + u_stretch
                if u_end <= n:
                    u_seq = seq[u_start:u_end]
                    u_count = u_seq.count('T')
                    
                    if u_count >= u_stretch * 0.75:
                        best_stem = matches
                        best_terminator = {
                            'start': left_start,
                            'end': u_end,
                            'stem_left': left,
                            'stem_right': right,
                            'loop': seq[left_end:right_start],
                            'u_stretch': u_seq,
                            'stem_matches': matches,
                            'stem_length': stem_length,
                            'loop_size': loop_size,
                            'u_count': u_count,
                            'score': matches + u_count,
                        }
        
        if best_terminator:
            terminators.append(best_terminator)
    
    terminators.sort(key=lambda x: x['score'], reverse=True)
    return terminators

def predict_rbs(seq, consensus='AGGAGG', min_score=3):
    if hasattr(seq, 'seq'):
        seq = seq.seq
    seq = seq.upper()
    
    rbs_sites = []
    n = len(seq)
    motif_len = len(consensus)
    
    for i in range(n - motif_len + 1):
        window = seq[i:i+motif_len]
        score = sum(1 for a, b in zip(window, consensus) if a == b)
        
        if score >= min_score:
            rbs_sites.append({
                'start': i,
                'end': i + motif_len,
                'sequence': window,
                'score': score,
                'consensus': consensus,
                'strand': '+',
            })
    
    rc = reverse_complement(seq)
    for i in range(n - motif_len + 1):
        window = rc[i:i+motif_len]
        score = sum(1 for a, b in zip(window, consensus) if a == b)
        
        if score >= min_score:
            orig_start = n - i - motif_len
            rbs_sites.append({
                'start': orig_start,
                'end': orig_start + motif_len,
                'sequence': window,
                'score': score,
                'consensus': consensus,
                'strand': '-',
            })
    
    rbs_sites.sort(key=lambda x: x['score'], reverse=True)
    return rbs_sites
