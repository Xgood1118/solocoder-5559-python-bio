import re

class Sequence:
    def __init__(self, seq, id='', description='', quality=None):
        self.seq = seq.upper() if isinstance(seq, str) else seq
        self.id = id
        self.description = description
        self.quality = quality
        self.annotations = {}
        self.features = []

    def __len__(self):
        return len(self.seq)

    def __str__(self):
        return self.seq

    def __getitem__(self, index):
        if isinstance(index, slice):
            return Sequence(self.seq[index], id=self.id, description=self.description)
        return self.seq[index]

    def __repr__(self):
        return f"Sequence(id='{self.id}', len={len(self.seq)})"

    @property
    def length(self):
        return len(self.seq)

def _parse_fasta_header(header):
    header = header.lstrip('>').strip()
    parts = header.split(None, 1)
    seq_id = parts[0] if parts else ''
    desc = parts[1] if len(parts) > 1 else ''
    return seq_id, desc

def read_fasta(file_path):
    sequences = []
    current_id = ''
    current_desc = ''
    current_seq = []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_seq:
                    sequences.append(Sequence(''.join(current_seq), current_id, current_desc))
                    current_seq = []
                current_id, current_desc = _parse_fasta_header(line)
            else:
                current_seq.append(line.upper())
        if current_seq:
            sequences.append(Sequence(''.join(current_seq), current_id, current_desc))
    return sequences

def read_fastq(file_path):
    sequences = []
    with open(file_path, 'r') as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq_line = f.readline().strip()
            plus_line = f.readline().strip()
            qual_line = f.readline().strip()
            if not header.startswith('@'):
                continue
            seq_id, desc = _parse_fasta_header(header.replace('@', '>', 1))
            quality = [ord(c) - 33 for c in qual_line]
            sequences.append(Sequence(seq_line.upper(), seq_id, desc, quality=quality))
    return sequences

def read_genbank(file_path):
    sequences = []
    with open(file_path, 'r') as f:
        content = f.read()
    records = re.split(r'^LOCUS\s+', content, flags=re.MULTILINE)
    for record in records[1:]:
        lines = ('LOCUS       ' + record).split('\n')
        locus_line = lines[0]
        locus_parts = locus_line.split()
        seq_id = locus_parts[1] if len(locus_parts) > 1 else ''
        definition = ''
        origin_start = -1
        features_start = -1
        for i, line in enumerate(lines):
            if line.startswith('DEFINITION'):
                definition = line[12:].strip()
                j = i + 1
                while j < len(lines) and lines[j].startswith(' ' * 12):
                    definition += ' ' + lines[j].strip()
                    j += 1
            if line.startswith('FEATURES'):
                features_start = i
            if line.startswith('ORIGIN'):
                origin_start = i
                break
        seq = ''
        if origin_start > 0:
            for line in lines[origin_start + 1:]:
                if line.startswith('//'):
                    break
                seq += re.sub(r'[\s\d/]', '', line).upper()
        features = []
        if features_start > 0 and origin_start > features_start:
            current_feature = None
            for line in lines[features_start + 1:origin_start]:
                if line.startswith(' ' * 5) and not line.startswith(' ' * 21):
                    if current_feature:
                        features.append(current_feature)
                    feat_match = re.match(r'\s{5}(\w+)\s+(.+)', line)
                    if feat_match:
                        current_feature = {
                            'type': feat_match.group(1),
                            'location': feat_match.group(2).strip(),
                            'qualifiers': {}
                        }
                elif current_feature and line.startswith(' ' * 21):
                    qual_match = re.match(r'\s{21}/(\w+)="?([^"]*)"?', line)
                    if qual_match:
                        key, value = qual_match.group(1), qual_match.group(2)
                        current_feature['qualifiers'][key] = value
            if current_feature:
                features.append(current_feature)
        seq_obj = Sequence(seq, id=seq_id, description=definition)
        seq_obj.features = features
        seq_obj.annotations['definition'] = definition
        sequences.append(seq_obj)
    return sequences

def read_clustal(file_path):
    sequences = []
    seq_dict = {}
    with open(file_path, 'r') as f:
        lines = f.readlines()
    for line in lines[1:]:
        line = line.strip()
        if not line or line.startswith(' '):
            continue
        if set(line) <= set('*:. '):
            continue
        parts = line.split()
        if len(parts) >= 2:
            name = parts[0]
            seq_part = parts[1]
            if name in seq_dict:
                seq_dict[name] += seq_part.upper()
            else:
                seq_dict[name] = seq_part.upper()
    for name, seq in seq_dict.items():
        sequences.append(Sequence(seq, id=name))
    return sequences

def write_fasta(sequences, file_path, line_width=60):
    with open(file_path, 'w') as f:
        for seq in sequences:
            header = seq.id
            if seq.description:
                header += ' ' + seq.description
            f.write(f'>{header}\n')
            for i in range(0, len(seq.seq), line_width):
                f.write(seq.seq[i:i+line_width] + '\n')

def write_genbank(sequences, file_path):
    with open(file_path, 'w') as f:
        for seq in sequences:
            locus = seq.id or 'UNKNOWN'
            f.write(f'LOCUS       {locus:<20} {len(seq.seq)} bp    DNA     linear   UNK\n')
            f.write(f'DEFINITION  {seq.description or "Unknown sequence."}\n')
            f.write(f'ORIGIN      \n')
            for i in range(0, len(seq.seq), 60):
                chunk = seq.seq[i:i+60].lower()
                pos = str(i + 1).rjust(9)
                formatted = ' '.join([chunk[j:j+10] for j in range(0, len(chunk), 10)])
                f.write(f'{pos} {formatted}\n')
            f.write('//\n')
