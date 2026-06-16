import re

class MiniBLAST:
    def __init__(self, database=None, word_size=11, score_threshold=30):
        self.database = {}
        self.word_size = word_size
        self.score_threshold = score_threshold
        if database:
            self.build_index(database)
    
    def build_index(self, sequences):
        self.database = {}
        self.seq_data = []
        
        for i, seq in enumerate(sequences):
            seq_str = seq.seq.upper() if hasattr(seq, 'seq') else seq.upper()
            name = seq.id if hasattr(seq, 'id') else f'seq_{i}'
            self.seq_data.append({'seq': seq_str, 'name': name, 'index': i})
            
            for j in range(len(seq_str) - self.word_size + 1):
                word = seq_str[j:j+self.word_size]
                if word not in self.database:
                    self.database[word] = []
                self.database[word].append((i, j))
    
    def search(self, query, top_n=10):
        query = query.upper()
        hits = []
        
        for i in range(len(query) - self.word_size + 1):
            word = query[i:i+self.word_size]
            if word in self.database:
                for (seq_idx, pos) in self.database[word]:
                    hits.append({
                        'seq_idx': seq_idx,
                        'query_start': i,
                        'subject_start': pos,
                        'word': word,
                    })
        
        if not hits:
            return []
        
        extended_hits = []
        for hit in hits:
            seq_idx = hit['seq_idx']
            subject_seq = self.seq_data[seq_idx]['seq']
            
            q_start = hit['query_start']
            s_start = hit['subject_start']
            
            q_left = q_start
            s_left = s_start
            left_score = 0
            left_matches = 0
            
            while q_left > 0 and s_left > 0 and left_score > -10:
                q_left -= 1
                s_left -= 1
                if query[q_left] == subject_seq[s_left]:
                    left_score += 2
                    left_matches += 1
                else:
                    left_score -= 1
            
            q_right = q_start + self.word_size
            s_right = s_start + self.word_size
            right_score = 0
            right_matches = 0
            
            while q_right < len(query) and s_right < len(subject_seq) and right_score > -10:
                if query[q_right] == subject_seq[s_right]:
                    right_score += 2
                    right_matches += 1
                else:
                    right_score -= 1
                q_right += 1
                s_right += 1
            
            total_score = self.word_size * 2 + left_score + right_score
            total_matches = self.word_size + left_matches + right_matches
            align_len = (q_right - q_left)
            
            if total_score >= self.score_threshold:
                extended_hits.append({
                    'subject': self.seq_data[seq_idx]['name'],
                    'subject_index': seq_idx,
                    'score': total_score,
                    'query_start': q_left,
                    'query_end': q_right,
                    'subject_start': s_left,
                    'subject_end': s_right,
                    'identity': total_matches / max(align_len, 1) * 100,
                    'align_length': align_len,
                    'query_align': query[q_left:q_right],
                    'subject_align': subject_seq[s_left:s_right],
                })
        
        if not extended_hits:
            return []
        
        extended_hits.sort(key=lambda x: x['score'], reverse=True)
        
        seen = set()
        unique_hits = []
        for hit in extended_hits:
            key = (hit['subject_index'], hit['query_start'], hit['subject_start'])
            if key not in seen:
                seen.add(key)
                unique_hits.append(hit)
                if len(unique_hits) >= top_n:
                    break
        
        return unique_hits
    
    def add_sequence(self, seq, name=None):
        seq_str = seq.seq.upper() if hasattr(seq, 'seq') else seq.upper()
        if name is None:
            name = f'seq_{len(self.seq_data)}'
        
        idx = len(self.seq_data)
        self.seq_data.append({'seq': seq_str, 'name': name, 'index': idx})
        
        for j in range(len(seq_str) - self.word_size + 1):
            word = seq_str[j:j+self.word_size]
            if word not in self.database:
                self.database[word] = []
            self.database[word].append((idx, j))
    
    def __len__(self):
        return len(self.seq_data)
