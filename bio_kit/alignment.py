import numpy as np

class AlignmentResult:
    def __init__(self, seq1_aligned, seq2_aligned, score, identity=0, gaps=0, align_len=0):
        self.seq1_aligned = seq1_aligned
        self.seq2_aligned = seq2_aligned
        self.score = score
        self.identity = identity
        self.gaps = gaps
        self.align_len = align_len

    def __str__(self):
        return f"Alignment(score={self.score}, identity={self.identity:.2f}%, gaps={self.gaps})"

    def get_match_line(self):
        match = ''
        for a, b in zip(self.seq1_aligned, self.seq2_aligned):
            if a == b and a != '-':
                match += '|'
            elif a != '-' and b != '-':
                match += '.'
            else:
                match += ' '
        return match

def _make_score_matrix(match=1, mismatch=-1, alphabet='ATCG'):
    matrix = {}
    for a in alphabet:
        matrix[a] = {}
        for b in alphabet:
            matrix[a][b] = match if a == b else mismatch
    return matrix

def needleman_wunsch(seq1, seq2, match_score=2, mismatch_penalty=-2, gap_penalty=-2):
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    
    n, m = len(seq1), len(seq2)
    dp = np.zeros((n + 1, m + 1), dtype=int)
    
    for i in range(n + 1):
        dp[i][0] = i * gap_penalty
    for j in range(m + 1):
        dp[0][j] = j * gap_penalty
    
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = dp[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty)
            delete = dp[i-1][j] + gap_penalty
            insert = dp[i][j-1] + gap_penalty
            dp[i][j] = max(match, delete, insert)
    
    align1 = ''
    align2 = ''
    i, j = n, m
    
    while i > 0 or j > 0:
        if i > 0 and j > 0 and dp[i][j] == dp[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty):
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1
            j -= 1
        elif i > 0 and dp[i][j] == dp[i-1][j] + gap_penalty:
            align1 = seq1[i-1] + align1
            align2 = '-' + align2
            i -= 1
        else:
            align1 = '-' + align1
            align2 = seq2[j-1] + align2
            j -= 1
    
    identity = sum(1 for a, b in zip(align1, align2) if a == b and a != '-')
    gaps = sum(1 for a in align1 if a == '-') + sum(1 for b in align2 if b == '-')
    align_len = len(align1)
    
    return AlignmentResult(
        align1, align2,
        score=dp[n][m],
        identity=identity / max(align_len, 1) * 100,
        gaps=gaps,
        align_len=align_len,
    )

def smith_waterman(seq1, seq2, match_score=2, mismatch_penalty=-2, gap_penalty=-2):
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    
    n, m = len(seq1), len(seq2)
    dp = np.zeros((n + 1, m + 1), dtype=int)
    
    max_score = 0
    max_i, max_j = 0, 0
    
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = dp[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty)
            delete = dp[i-1][j] + gap_penalty
            insert = dp[i][j-1] + gap_penalty
            dp[i][j] = max(0, match, delete, insert)
            if dp[i][j] > max_score:
                max_score = dp[i][j]
                max_i, max_j = i, j
    
    align1 = ''
    align2 = ''
    i, j = max_i, max_j
    
    while i > 0 and j > 0 and dp[i][j] > 0:
        if dp[i][j] == dp[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty):
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1
            j -= 1
        elif dp[i][j] == dp[i-1][j] + gap_penalty:
            align1 = seq1[i-1] + align1
            align2 = '-' + align2
            i -= 1
        else:
            align1 = '-' + align1
            align2 = seq2[j-1] + align2
            j -= 1
    
    identity = sum(1 for a, b in zip(align1, align2) if a == b and a != '-')
    gaps = sum(1 for a in align1 if a == '-') + sum(1 for b in align2 if b == '-')
    align_len = len(align1)
    
    result = AlignmentResult(
        align1, align2,
        score=max_score,
        identity=identity / max(align_len, 1) * 100,
        gaps=gaps,
        align_len=align_len,
    )
    result.start1 = i
    result.start2 = j
    result.end1 = max_i
    result.end2 = max_j
    
    return result

def _pairwise_distance(seqs):
    n = len(seqs)
    dist = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            aln = needleman_wunsch(seqs[i], seqs[j])
            identity = aln.identity / 100.0
            d = 1.0 - identity
            dist[i][j] = d
            dist[j][i] = d
    return dist

def _profile_consensus(profile):
    consensus = []
    for col in profile:
        max_count = max(col.values())
        for nuc in 'ATCG':
            if col.get(nuc, 0) == max_count:
                consensus.append(nuc)
                break
    return ''.join(consensus)

def _align_to_profile(seq, profile):
    consensus = _profile_consensus(profile)
    aln = needleman_wunsch(seq, consensus)
    
    new_profile = []
    gap_count = profile[0].get('-', 0) if profile else 0
    total = sum(profile[0].values()) if profile else 0
    
    seq_pos = 0
    prof_pos = 0
    
    for a, b in zip(aln.seq1_aligned, aln.seq2_aligned):
        col = {}
        if a != '-':
            col[a] = col.get(a, 0) + 1
            seq_pos += 1
        if b != '-' and prof_pos < len(profile):
            for k, v in profile[prof_pos].items():
                col[k] = col.get(k, 0) + v
            prof_pos += 1
        if a == '-' and prof_pos < len(profile):
            for k, v in profile[prof_pos].items():
                col[k] = col.get(k, 0) + v
            col['-'] = col.get('-', 0) + 1
        elif b == '-':
            col['-'] = col.get('-', 0) + total
        
        new_profile.append(col)
    
    while prof_pos < len(profile):
        col = dict(profile[prof_pos])
        col['-'] = col.get('-', 0) + 1
        new_profile.append(col)
        prof_pos += 1
    
    return new_profile, aln.seq1_aligned

def _seq_to_profile(seq):
    profile = []
    for c in seq:
        profile.append({c: 1})
    return profile

def _merge_profiles(profile1, profile2):
    pass

def clustalw_style_msa(sequences, match_score=2, mismatch_penalty=-2, gap_penalty=-2):
    if isinstance(sequences[0], str):
        seq_list = [s.upper() for s in sequences]
        names = [f'seq{i+1}' for i in range(len(seq_list))]
    else:
        seq_list = [s.seq.upper() for s in sequences]
        names = [s.id for s in sequences]
    
    n = len(seq_list)
    if n < 2:
        return [seq_list[0]]
    
    dist = _pairwise_distance(seq_list)
    
    clusters = [[i] for i in range(n)]
    newick_parts = [names[i] for i in range(n)]
    
    while len(clusters) > 1:
        min_dist = float('inf')
        best_i, best_j = 0, 1
        
        for i in range(len(clusters)):
            for j in range(i + 1, len(clusters)):
                total = 0
                count = 0
                for ci in clusters[i]:
                    for cj in clusters[j]:
                        total += dist[ci][cj]
                        count += 1
                avg = total / count if count > 0 else 0
                if avg < min_dist:
                    min_dist = avg
                    best_i, best_j = i, j
        
        new_cluster = clusters[best_i] + clusters[best_j]
        
        left = newick_parts[best_i]
        right = newick_parts[best_j]
        newick_parts[best_i] = f'({left}:{min_dist/2},{right}:{min_dist/2})'
        newick_parts.pop(best_j)
        
        clusters[best_i] = new_cluster
        clusters.pop(best_j)
    
    final_newick = newick_parts[0] + ';'
    
    order = []
    def parse_newick(newick_str):
        newick_str = newick_str.rstrip(';')
        def parse_recursive(s, idx):
            if s[idx] != '(':
                name_end = idx
                while name_end < len(s) and s[name_end] not in '(),:':
                    name_end += 1
                name = s[idx:name_end]
                return name, name_end
            else:
                idx += 1
                children = []
                while True:
                    child, idx = parse_recursive(s, idx)
                    if idx < len(s) and s[idx] == ':':
                        idx += 1
                        while idx < len(s) and s[idx] not in '(),':
                            idx += 1
                    children.append(child)
                    if idx < len(s) and s[idx] == ',':
                        idx += 1
                    elif idx < len(s) and s[idx] == ')':
                        idx += 1
                        break
                return children, idx
        tree, _ = parse_recursive(newick_str, 0)
        def flatten(t):
            if isinstance(t, str):
                order.append(t)
            else:
                for c in t:
                    flatten(c)
        flatten(tree)
    
    parse_newick(final_newick)
    
    ordered_indices = [names.index(name) for name in order if name in names]
    
    aligned = {}
    if n >= 2:
        first = ordered_indices[0]
        second = ordered_indices[1]
        aln = needleman_wunsch(seq_list[first], seq_list[second])
        profile = _seq_to_profile(aln.seq1_aligned)
        aligned[names[first]] = aln.seq1_aligned
        aligned[names[second]] = aln.seq2_aligned
        
        for idx in ordered_indices[2:]:
            name = names[idx]
            seq = seq_list[idx]
            consensus = _profile_consensus(profile)
            aln = needleman_wunsch(seq, consensus)
            
            new_profile = []
            s_pos = 0
            c_pos = 0
            
            for a, b in zip(aln.seq1_aligned, aln.seq2_aligned):
                col = {}
                if a != '-':
                    col[a] = 1
                    s_pos += 1
                else:
                    col['-'] = 1
                
                if b != '-' and c_pos < len(profile):
                    for k, v in profile[c_pos].items():
                        col[k] = col.get(k, 0) + v
                    c_pos += 1
                elif b == '-':
                    total = sum(profile[0].values()) if profile else 1
                    col['-'] = col.get('-', 0) + total
                
                new_profile.append(col)
            
            while c_pos < len(profile):
                col = dict(profile[c_pos])
                col['-'] = col.get('-', 0) + 1
                new_profile.append(col)
                c_pos += 1
            
            profile = new_profile
            
            for old_name in list(aligned.keys()):
                old_seq = aligned[old_name]
                new_old = ''
                p = 0
                for b in aln.seq2_aligned:
                    if b != '-':
                        if p < len(old_seq):
                            new_old += old_seq[p]
                            p += 1
                    else:
                        new_old += '-'
                while p < len(old_seq):
                    new_old += old_seq[p]
                    p += 1
                aligned[old_name] = new_old
            
            aligned[name] = aln.seq1_aligned
    
    result = []
    for name in order:
        if name in aligned:
            result.append(aligned[name])
    
    class MSAResult:
        def __init__(self, aligned_seqs, seq_names, newick):
            self.aligned_sequences = aligned_seqs
            self.names = seq_names
            self.newick = newick
            self.num_seqs = len(aligned_seqs)
            self.alignment_length = len(aligned_seqs[0]) if aligned_seqs else 0
        
        def __repr__(self):
            return f"MSAResult(num_seqs={self.num_seqs}, length={self.alignment_length})"
        
        def get_consensus(self, threshold=0.5):
            if not self.aligned_sequences:
                return ''
            consensus = []
            n = self.num_seqs
            for i in range(self.alignment_length):
                counts = {}
                for seq in self.aligned_sequences:
                    if i < len(seq):
                        c = seq[i]
                        counts[c] = counts.get(c, 0) + 1
                max_count = max(counts.values())
                best_char = max(counts, key=counts.get)
                if max_count / n >= threshold:
                    consensus.append(best_char)
                else:
                    consensus.append('X' if best_char != '-' else '-')
            return ''.join(consensus)
        
        def identity_matrix(self):
            n = self.num_seqs
            matrix = [[0.0] * n for _ in range(n)]
            for i in range(n):
                matrix[i][i] = 100.0
                for j in range(i + 1, n):
                    same = sum(1 for a, b in zip(self.aligned_sequences[i], self.aligned_sequences[j]) 
                              if a == b and a != '-')
                    total = sum(1 for a, b in zip(self.aligned_sequences[i], self.aligned_sequences[j]) 
                               if a != '-' or b != '-')
                    matrix[i][j] = same / total * 100 if total > 0 else 0
                    matrix[j][i] = matrix[i][j]
            return matrix
    
    return MSAResult(result, order, final_newick)
