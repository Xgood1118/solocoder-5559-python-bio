import math
import io

def _gc_color(gc):
    if gc < 30:
        r = int(0 + gc * 3.4)
        g = int(100 + gc * 4.3)
        b = int(255)
    elif gc < 50:
        t = (gc - 30) / 20
        r = int(100 + t * 100)
        g = int(186 + t * 50)
        b = int(255 - t * 128)
    else:
        t = (gc - 50) / 50
        r = int(200 + t * 55)
        g = int(236 - t * 100)
        b = int(127 - t * 127)
    return f'rgb({r}, {g}, {b})'

def draw_sequence_gc(seq, window=50, width=800, height=200, output_file=None):
    if hasattr(seq, 'seq'):
        seq = seq.seq
    seq = seq.upper()
    n = len(seq)
    
    gc_values = []
    for i in range(0, n - window + 1, max(1, window // 10)):
        window_seq = seq[i:i+window]
        gc = (window_seq.count('G') + window_seq.count('C')) / window * 100
        gc_values.append((i, gc))
    
    margin = 60
    plot_width = width - 2 * margin
    plot_height = height - 2 * margin
    
    max_gc = 100
    min_gc = 0
    
    points = []
    for pos, gc in gc_values:
        x = margin + pos / max(n - window, 1) * plot_width
        y = margin + plot_height - (gc - min_gc) / (max_gc - min_gc) * plot_height
        points.append((x, y))
    
    path_d = f'M {points[0][0]:.1f} {points[0][1]:.1f}'
    for x, y in points[1:]:
        path_d += f' L {x:.1f} {y:.1f}'
    
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '  <rect width="100%" height="100%" fill="white"/>',
        '  <text x="10" y="20" font-size="14" font-family="Arial" fill="black">GC Content (%)</text>',
    ]
    
    for gc in range(0, 101, 25):
        y = margin + plot_height - gc / (max_gc - min_gc) * plot_height
        svg.append(f'  <line x1="{margin}" y1="{y}" x2="{width - margin}" y2="{y}" stroke="#ddd" stroke-width="1"/>')
        svg.append(f'  <text x="{margin - 5}" y="{y + 4}" font-size="10" text-anchor="end" fill="#666">{gc}</text>')
    
    svg.append(f'  <path d="{path_d}" fill="none" stroke="#2c7fb8" stroke-width="2"/>')
    fill_path = path_d + f' L {points[-1][0]:.1f} {margin + plot_height} L {points[0][0]:.1f} {margin + plot_height} Z'
    svg.append(f'  <path d="{fill_path}" fill="#2c7fb8" fill-opacity="0.2" stroke="none"/>')
    
    svg.append(f'  <text x="{width // 2}" y="{height - 20}" font-size="12" text-anchor="middle" fill="#666">Position (bp)</text>')
    svg.append(f'  <text x="{margin}" y="{height - 5}" font-size="10" fill="#666">0</text>')
    svg.append(f'  <text x="{width - margin}" y="{height - 5}" font-size="10" text-anchor="end" fill="#666">{n}</text>')
    
    svg.append('</svg>')
    
    svg_content = '\n'.join(svg)
    
    if output_file:
        with open(output_file, 'w') as f:
            f.write(svg_content)
    
    return svg_content

def gc_plot(seq, window=50, width=800, height=200, output_file=None):
    return draw_sequence_gc(seq, window, width, height, output_file)

def alignment_colors(aligned_sequences, names=None, width=800, char_width=8, line_height=20, output_file=None):
    if names is None:
        names = [f'seq{i+1}' for i in range(len(aligned_sequences))]
    
    n_seqs = len(aligned_sequences)
    if n_seqs == 0:
        return ''
    
    align_len = len(aligned_sequences[0])
    
    color_map = {
        'A': '#88e088',
        'T': '#88b8f0',
        'G': '#f0e080',
        'C': '#f08888',
        'U': '#88b8f0',
        '-': '#f0f0f0',
        'N': '#cccccc',
    }
    
    protein_colors = {
        'A': '#88e088', 'V': '#88e088', 'L': '#88e088', 'I': '#88e088', 'M': '#88e088',
        'K': '#88b8f0', 'R': '#88b8f0', 'H': '#88b8f0',
        'D': '#f08888', 'E': '#f08888',
        'S': '#f0e080', 'T': '#f0e080',
        'N': '#c0a0e0', 'Q': '#c0a0e0',
        'C': '#f0a040',
        'W': '#e080e0',
        'Y': '#e080e0',
        'F': '#e080e0',
        'P': '#a0a0a0',
        'G': '#d0d0d0',
        '*': '#ff0000',
        '-': '#f0f0f0',
        'X': '#cccccc',
    }
    
    sample = aligned_sequences[0]
    is_protein = any(c.upper() in 'VLIKRHDESTNQCWYFPGX' for c in sample)
    
    colors = protein_colors if is_protein else color_map
    
    margin_left = 120
    margin_top = 30
    total_width = margin_left + align_len * char_width
    total_height = margin_top + n_seqs * line_height + 30
    
    chars_per_line = min(align_len, (width - margin_left) // char_width)
    total_lines = math.ceil(align_len / chars_per_line) if chars_per_line > 0 else 1
    total_height = margin_top + n_seqs * line_height * total_lines + 30 * total_lines
    
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{total_height}" viewBox="0 0 {width} {total_height}" font-family="monospace">',
        '  <rect width="100%" height="100%" fill="white"/>',
        f'  <text x="10" y="15" font-size="12" fill="black" font-weight="bold">Sequence Alignment</text>',
    ]
    
    for line_idx in range(total_lines):
        start_col = line_idx * chars_per_line
        end_col = min(start_col + chars_per_line, align_len)
        
        y_offset = margin_top + line_idx * (n_seqs * line_height + 30)
        
        svg.append(f'  <text x="{margin_left}" y="{y_offset - 5}" font-size="10" fill="#666">{start_col + 1}</text>')
        svg.append(f'  <text x="{margin_left + (end_col - start_col) * char_width}" y="{y_offset - 5}" font-size="10" text-anchor="end" fill="#666">{end_col}</text>')
        
        for seq_idx in range(n_seqs):
            y = y_offset + seq_idx * line_height
            name = names[seq_idx]
            display_name = name[:14] + '..' if len(name) > 14 else name
            svg.append(f'  <text x="5" y="{y + 14}" font-size="11" fill="#333" font-family="Arial">{display_name}</text>')
            
            for col in range(start_col, end_col):
                x = margin_left + (col - start_col) * char_width
                char = aligned_sequences[seq_idx][col].upper()
                bg = colors.get(char, '#ffffff')
                
                svg.append(f'  <rect x="{x}" y="{y}" width="{char_width}" height="{line_height - 2}" fill="{bg}"/>')
                svg.append(f'  <text x="{x + char_width // 2}" y="{y + 14}" font-size="10" text-anchor="middle" fill="#333">{char}</text>')
    
    svg.append('</svg>')
    
    svg_content = '\n'.join(svg)
    
    if output_file:
        with open(output_file, 'w') as f:
            f.write(svg_content)
    
    return svg_content

def tree_to_svg(tree, width=600, height=400, output_file=None):
    import os
    import tempfile
    import sys

    ete_tree = tree.to_ete3_tree() if hasattr(tree, 'to_ete3_tree') else None

    ete_bootstrap_info = {}
    if ete_tree is not None:
        try:
            leaf_index = 0
            for node in ete_tree.traverse("preorder"):
                if node.is_leaf():
                    leaf_index += 1
                elif node.support and node.support > 0:
                    child_leaves = tuple(sorted(node.get_leaf_names()))
                    ete_bootstrap_info[child_leaves] = float(node.support)
        except Exception:
            pass

    if ete_tree is not None and sys.platform != 'win32':
        try:
            os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
            tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.svg', delete=False)
            tmp.close()
            try:
                from ete3 import NodeStyle, TextFace, AttrFace

                for node in ete_tree.traverse():
                    ns = NodeStyle()
                    ns["size"] = 0
                    ns["vt_line_width"] = 2
                    ns["hz_line_width"] = 2
                    node.set_style(ns)

                    if not node.is_leaf() and node.support and node.support > 0:
                        support_face = TextFace(" %d" % int(node.support), fsize=11, fgcolor="#c0392b")
                        node.add_face(support_face, column=0, position="branch-top")

                ete_tree.render(tmp.name, w=width, h=height, units='px')
                with open(tmp.name, 'r', encoding='utf-8', errors='ignore') as f:
                    svg_content = f.read()
                if output_file:
                    with open(output_file, 'w', encoding='utf-8') as f:
                        f.write(svg_content)
                return svg_content
            finally:
                try:
                    os.unlink(tmp.name)
                except Exception:
                    pass
        except Exception:
            pass

    leaves = tree.get_leaves()
    n_leaves = len(leaves)

    leaf_y_step = (height - 60) / max(n_leaves, 1)
    margin_left = 80
    margin_right = 30
    margin_top = 30
    margin_bottom = 30

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" font-family="Arial">',
        '  <rect width="100%" height="100%" fill="white"/>',
        f'  <text x="{width // 2}" y="20" font-size="14" text-anchor="middle" font-weight="bold">Phylogenetic Tree</text>',
    ]

    def calc_total_distance(node):
        if node.is_leaf():
            return node.distance
        max_dist = 0
        for child in node.children:
            d = calc_total_distance(child)
            max_dist = max(max_dist, d + node.distance)
        return max_dist

    total_dist = 0
    def calc_max_depth(node, depth=0):
        nonlocal total_dist
        if node.is_leaf():
            total_dist = max(total_dist, depth + node.distance)
            return
        for child in node.children:
            calc_max_depth(child, depth + node.distance)

    if tree.root:
        calc_max_depth(tree.root)

    if total_dist == 0:
        total_dist = 1.0

    leaf_index = [0]

    def draw_node(node, x, y, depth):
        if node.is_leaf():
            y_pos = margin_top + leaf_index[0] * leaf_y_step + leaf_y_step / 2
            leaf_index[0] += 1

            x_end = margin_left + depth * ((width - margin_left - margin_right) / total_dist)
            svg.append(f'  <line x1="{x}" y1="{y}" x2="{x_end}" y2="{y_pos}" stroke="#333" stroke-width="1.5"/>')
            svg.append(f'  <text x="{x_end + 5}" y="{y_pos + 4}" font-size="11" fill="#333">{node.name}</text>')

            return y_pos, x_end, {node.name}
        else:
            child_ys = []
            child_xs = []
            all_leaves = set()

            new_depth = depth + node.distance
            x_pos = margin_left + new_depth * ((width - margin_left - margin_right) / total_dist)

            for child in node.children:
                cy, cx, cleaves = draw_node(child, x_pos, 0, new_depth)
                child_ys.append(cy)
                child_xs.append(cx)
                all_leaves |= cleaves

            avg_y = sum(child_ys) / len(child_ys)

            svg.append(f'  <line x1="{x}" y1="{avg_y}" x2="{x_pos}" y2="{avg_y}" stroke="#333" stroke-width="1.5"/>')

            if len(child_ys) >= 2:
                min_y = min(child_ys)
                max_y = max(child_ys)
                svg.append(f'  <line x1="{x_pos}" y1="{min_y}" x2="{x_pos}" y2="{max_y}" stroke="#333" stroke-width="1.5"/>')

            bootstrap_val = node.bootstrap
            if bootstrap_val is not None:
                leaves_key = tuple(sorted(all_leaves))
                if leaves_key in ete_bootstrap_info:
                    bootstrap_val = ete_bootstrap_info[leaves_key]

            if bootstrap_val is not None and not node.is_leaf():
                svg.append(f'  <circle cx="{x_pos}" cy="{avg_y}" r="5" fill="white" stroke="#c0392b" stroke-width="1.5"/>')
                svg.append(f'  <text x="{x_pos + 8}" y="{avg_y - 8}" font-size="11" fill="#c0392b" font-weight="bold">{bootstrap_val:.0f}</text>')

            return avg_y, x_pos, all_leaves

    if tree.root:
        draw_node(tree.root, margin_left, 0, 0)

    svg.append('</svg>')

    svg_content = '\n'.join(svg)

    if output_file:
        with open(output_file, 'w') as f:
            f.write(svg_content)

    return svg_content

def similarity_heatmap(sequences, names=None, width=500, height=400, output_file=None):
    from .alignment import needleman_wunsch
    
    if isinstance(sequences[0], str):
        seqs = [s.upper() for s in sequences]
    else:
        seqs = [s.seq.upper() for s in sequences]
    
    if names is None:
        names = [f'seq{i+1}' for i in range(len(seqs))]
    
    n = len(seqs)
    
    matrix = []
    for i in range(n):
        row = []
        for j in range(n):
            if i == j:
                row.append(100.0)
            else:
                aln = needleman_wunsch(seqs[i], seqs[j])
                row.append(aln.identity)
        matrix.append(row)
    
    margin_left = 100
    margin_top = 40
    margin_right = 30
    margin_bottom = 100
    
    cell_width = (width - margin_left - margin_right) / max(n, 1)
    cell_height = (height - margin_top - margin_bottom) / max(n, 1)
    
    def heat_color(value):
        v = value / 100.0
        if v < 0.5:
            r = int(255 * (1 - v * 2))
            g = int(255 * v * 2)
            b = 255
        else:
            r = 0
            g = 255
            b = int(255 * (1 - v) * 2)
        return f'rgb({r}, {g}, {b})'
    
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" font-family="Arial">',
        '  <rect width="100%" height="100%" fill="white"/>',
        f'  <text x="{width // 2}" y="20" font-size="14" text-anchor="middle" font-weight="bold">Sequence Similarity Matrix (%)</text>',
    ]
    
    for i in range(n):
        for j in range(n):
            x = margin_left + j * cell_width
            y = margin_top + i * cell_height
            value = matrix[i][j]
            color = heat_color(value)
            
            svg.append(f'  <rect x="{x}" y="{y}" width="{cell_width:.1f}" height="{cell_height:.1f}" fill="{color}" stroke="white" stroke-width="0.5"/>')
            
            if cell_width > 30 and cell_height > 20:
                svg.append(f'  <text x="{x + cell_width / 2}" y="{y + cell_height / 2 + 4}" font-size="10" text-anchor="middle" fill="white" font-weight="bold">{value:.0f}</text>')
    
    for i in range(n):
        name = names[i][:15] + '..' if len(names[i]) > 15 else names[i]
        y = margin_top + i * cell_height + cell_height / 2 + 4
        svg.append(f'  <text x="{margin_left - 5}" y="{y}" font-size="10" text-anchor="end" fill="#333">{name}</text>')
        
        x = margin_left + i * cell_width + cell_width / 2
        svg.append(f'  <text x="{x}" y="{height - margin_bottom + 15}" font-size="10" text-anchor="end" fill="#333" transform="rotate(-45, {x}, {height - margin_bottom + 15})">{name}</text>')
    
    svg.append('</svg>')
    
    svg_content = '\n'.join(svg)
    
    if output_file:
        with open(output_file, 'w') as f:
            f.write(svg_content)
    
    return svg_content
