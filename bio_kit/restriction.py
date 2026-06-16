import re
from .data.enzymes import RESTRICTION_ENZYMES

class RestrictionEnzyme:
    def __init__(self, name, recognition, cut_top, cut_bottom, enzyme_type='II'):
        self.name = name
        self.recognition = recognition.upper()
        self.cut_top = cut_top
        self.cut_bottom = cut_bottom
        self.type = enzyme_type

    def __repr__(self):
        return f"RestrictionEnzyme('{self.name}', '{self.recognition}')"

    def __str__(self):
        return f"{self.name}: {self.recognition}"

def _pattern_to_regex(pattern):
    iupac = {
        'A': 'A', 'T': 'T', 'G': 'G', 'C': 'C',
        'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
        'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
        'H': '[ACT]', 'V': '[ACG]', 'N': '[ATCG]',
    }
    regex = ''
    i = 0
    while i < len(pattern):
        c = pattern[i]
        if c == '(':
            j = i + 1
            while j < len(pattern) and pattern[j] != ')':
                j += 1
            nums = pattern[i+1:j].split(',')
            if len(nums) == 1:
                regex += '{' + nums[0] + '}'
            else:
                regex += '{' + nums[0] + ',' + nums[1] + '}'
            i = j + 1
        else:
            regex += iupac.get(c, c)
            i += 1
    return regex

def _count_pattern_len(pattern):
    i = 0
    length = 0
    while i < len(pattern):
        c = pattern[i]
        if c == '(':
            j = i + 1
            while j < len(pattern) and pattern[j] != ')':
                j += 1
            nums = pattern[i+1:j].split(',')
            length += int(nums[0])
            i = j + 1
        elif c in 'ATGCRYSWKMBDHVN':
            length += 1
            i += 1
        else:
            i += 1
    return length

def get_enzyme(name):
    info = RESTRICTION_ENZYMES.get(name)
    if info:
        return RestrictionEnzyme(name, info['recognition'], info['cut_top'], info['cut_bottom'], info.get('type', 'II'))
    return None

def list_enzymes():
    return list(RESTRICTION_ENZYMES.keys())

def find_restriction_sites(seq, enzymes=None):
    if enzymes is None:
        enzymes = list(RESTRICTION_ENZYMES.keys())
    if isinstance(seq, str):
        seq = seq.upper()
    else:
        seq = seq.seq.upper()
    
    results = {}
    for enzyme_name in enzymes:
        info = RESTRICTION_ENZYMES.get(enzyme_name)
        if not info:
            continue
        pattern = info['recognition']
        cut_top = info['cut_top']
        cut_bottom = info['cut_bottom']
        
        regex_pattern = _pattern_to_regex(pattern)
        sites = []
        for match in re.finditer(r'(?=(' + regex_pattern + '))', seq):
            pos = match.start()
            sites.append({
                'position': pos,
                'cut_top': pos + cut_top,
                'cut_bottom': pos + cut_bottom,
                'match': match.group(1),
            })
        
        if sites:
            fragment_lengths = []
            cuts = sorted([s['cut_top'] for s in sites])
            prev = 0
            for cut in cuts:
                fragment_lengths.append(cut - prev)
                prev = cut
            fragment_lengths.append(len(seq) - prev)
            
            results[enzyme_name] = {
                'sites': sites,
                'count': len(sites),
                'fragments': fragment_lengths,
            }
    
    return results
