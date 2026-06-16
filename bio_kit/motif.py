import re
from .data.prosite import PROSITE_PATTERNS, prosite_to_regex

def find_motif_regex(seq, pattern):
    if hasattr(seq, 'seq'):
        seq = seq.seq
    seq = seq.upper()
    
    matches = []
    for match in re.finditer(pattern, seq):
        matches.append({
            'start': match.start(),
            'end': match.end(),
            'match': match.group(),
            'pattern': pattern,
        })
    return matches

def find_prosite_motif(seq, prosite_id=None, pattern=None):
    if hasattr(seq, 'seq'):
        seq = seq.seq
    seq = seq.upper()
    
    if pattern is None and prosite_id is not None:
        info = PROSITE_PATTERNS.get(prosite_id)
        if info:
            pattern = info['pattern']
        else:
            return []
    
    if pattern is None:
        return []
    
    regex_pattern = prosite_to_regex(pattern)
    
    matches = []
    for match in re.finditer(regex_pattern, seq):
        matches.append({
            'start': match.start(),
            'end': match.end(),
            'match': match.group(),
            'prosite_id': prosite_id,
            'pattern': pattern,
            'name': PROSITE_PATTERNS.get(prosite_id, {}).get('name', ''),
            'description': PROSITE_PATTERNS.get(prosite_id, {}).get('description', ''),
        })
    return matches

def find_all_prosite_motifs(seq):
    if hasattr(seq, 'seq'):
        seq = seq.seq
    seq = seq.upper()
    
    all_matches = []
    for ps_id, info in PROSITE_PATTERNS.items():
        try:
            matches = find_prosite_motif(seq, prosite_id=ps_id)
            all_matches.extend(matches)
        except:
            continue
    return all_matches
