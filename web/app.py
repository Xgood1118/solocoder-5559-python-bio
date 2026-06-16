import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from flask import Flask, render_template, request, jsonify, send_file, Response
import io
import tempfile

app = Flask(__name__)
app.config['SECRET_KEY'] = 'bio_kit_secret_key'

from bio_kit import (
    read_fasta, read_fastq, read_genbank, read_clustal,
    Sequence, reverse_complement, complement, transcribe, translate, find_orfs, six_frame_translation,
    find_restriction_sites, get_enzyme, list_enzymes,
    design_primers, check_dimer, check_hairpin, calculate_tm, calculate_gc,
    needleman_wunsch, smith_waterman, clustalw_style_msa,
    MiniBLAST,
    upgma, neighbor_joining, maximum_parsimony, bootstrap, DistanceModel,
    find_motif_regex, find_prosite_motif, PROSITE_PATTERNS,
    gc_content, n50, kmer_frequency, shannon_entropy, rscu_analysis, sequence_stats,
    predict_promoters, predict_terminators, predict_rbs,
    gc_plot, alignment_colors, tree_to_svg, similarity_heatmap,
    export_fasta, export_newick,
)

EXAMPLE_DNA = """
>Sample DNA sequence
ATGGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT
CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA
TCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC
GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT
CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA
TCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
TAG
""".strip()

EXAMPLE_MULTI = """
>Sequence_A
ATGGCGATTACCGTTGATGTTGATGCTGGTGAAACTGCTGGTGAAACCGGTGAA
>Sequence_B
ATGGCGATTACCGTTCATGTTGATGCTGGTCAAACTGCTGGTGAAACCGGTGAA
>Sequence_C
ATGGCGATTACCGTTGACGTTAATGCTGGTGAAACGGCTGGCGAAACCGGCAAA
>Sequence_D
ATGGCTATCACCGTTGATGTTAATGCAGGTGAAACTGCTGGCGAGACTGGTGAA
""".strip()


def _parse_input(text, fmt='fasta'):
    if not text.strip():
        return []
    
    tmp = tempfile.NamedTemporaryFile(mode='w', suffix=f'.{fmt}', delete=False)
    tmp.write(text)
    tmp.close()
    
    try:
        if fmt == 'fasta':
            seqs = read_fasta(tmp.name)
        elif fmt == 'fastq':
            seqs = read_fastq(tmp.name)
        elif fmt == 'genbank':
            seqs = read_genbank(tmp.name)
        elif fmt == 'clustal':
            seqs = read_clustal(tmp.name)
        else:
            seqs = read_fasta(tmp.name)
    finally:
        os.unlink(tmp.name)
    
    return seqs


@app.route('/')
def index():
    return render_template('index.html')


@app.route('/api/analyze', methods=['POST'])
def api_analyze():
    data = request.json
    seq_text = data.get('sequence', '')
    seq_format = data.get('format', 'fasta')
    analysis = data.get('analysis', 'basic')
    
    try:
        seqs = _parse_input(seq_text, seq_format)
    except Exception as e:
        return jsonify({'error': str(e)}), 400
    
    if not seqs:
        return jsonify({'error': 'No sequences found'}), 400
    
    result = {}
    
    if analysis == 'basic':
        seq = seqs[0]
        result['sequence'] = seq.seq
        result['length'] = len(seq)
        result['gc'] = gc_content(seq)
        result['reverse_complement'] = reverse_complement(seq)
        result['transcribed'] = transcribe(seq)
        result['translated'] = translate(seq)
        stats = sequence_stats(seq)
        result['stats'] = stats
        result['entropy'] = shannon_entropy(seq)
    
    elif analysis == 'orf':
        seq = seqs[0]
        orfs = find_orfs(seq, min_length=30)
        result['orfs'] = orfs[:20]
        result['total_orfs'] = len(orfs)
        result['six_frames'] = six_frame_translation(seq)
    
    elif analysis == 'restriction':
        seq = seqs[0]
        enzymes = data.get('enzymes', ['EcoRI', 'BamHI', 'HindIII'])
        sites = find_restriction_sites(seq.seq, enzymes)
        result['sites'] = sites
        result['available_enzymes'] = list_enzymes()[:50]
    
    elif analysis == 'primer':
        seq = seqs[0]
        start = int(data.get('start', 0))
        end = int(data.get('end', len(seq)))
        min_len = int(data.get('min_length', 18))
        max_len = int(data.get('max_length', 25))
        min_tm = float(data.get('min_tm', 55))
        max_tm = float(data.get('max_tm', 65))
        min_gc = float(data.get('min_gc', 40))
        max_gc = float(data.get('max_gc', 60))
        
        primers = design_primers(
            seq.seq, start=start, end=end,
            length_range=(min_len, max_len),
            tm_range=(min_tm, max_tm),
            gc_range=(min_gc, max_gc),
        )
        result['primers'] = primers
    
    elif analysis == 'alignment_pair':
        if len(seqs) < 2:
            return jsonify({'error': 'Need at least 2 sequences'}), 400
        seq1 = seqs[0].seq
        seq2 = seqs[1].seq
        algo = data.get('algorithm', 'needleman')
        
        if algo == 'needleman':
            aln = needleman_wunsch(seq1, seq2)
        else:
            aln = smith_waterman(seq1, seq2)
        
        result['alignment'] = {
            'seq1': aln.seq1_aligned,
            'seq2': aln.seq2_aligned,
            'match_line': aln.get_match_line(),
            'score': aln.score,
            'identity': aln.identity,
            'gaps': aln.gaps,
            'length': aln.align_len,
            'name1': seqs[0].id or 'seq1',
            'name2': seqs[1].id or 'seq2',
        }
    
    elif analysis == 'alignment_multi':
        if len(seqs) < 2:
            return jsonify({'error': 'Need at least 2 sequences'}), 400
        seqs_list = [s.seq for s in seqs]
        names = [s.id or f'seq{i+1}' for i, s in enumerate(seqs)]
        
        msa = clustalw_style_msa(seqs_list)
        result['msa'] = {
            'sequences': msa.aligned_sequences,
            'names': msa.names,
            'consensus': msa.get_consensus(),
            'length': msa.alignment_length,
            'num_seqs': msa.num_seqs,
            'newick': msa.newick,
        }
    
    elif analysis == 'blast':
        if len(seqs) < 2:
            return jsonify({'error': 'Need at least 2 sequences (query + database)'}), 400
        
        query = seqs[0].seq
        db_seqs = seqs[1:]
        
        blast = MiniBLAST(db_seqs, word_size=7, score_threshold=20)
        results = blast.search(query, top_n=10)
        result['blast'] = {
            'hits': results,
            'db_size': len(db_seqs),
            'query_length': len(query),
        }
    
    elif analysis == 'phylogeny':
        if len(seqs) < 3:
            return jsonify({'error': 'Need at least 3 sequences for phylogeny'}), 400
        
        seqs_list = [s.seq for s in seqs]
        names = [s.id or f'seq{i+1}' for i, s in enumerate(seqs)]
        method = data.get('method', 'upgma')
        dist_model = data.get('distance_model', 'jukes_cantor')
        do_bootstrap = data.get('bootstrap', False)
        n_boot = int(data.get('n_bootstrap', 100))
        
        if method == 'upgma':
            if do_bootstrap:
                tree = bootstrap(seqs_list, names, method='upgma', 
                               n_bootstraps=n_boot, distance_model=dist_model)
            else:
                tree = upgma(seqs_list, names, distance_model=dist_model)
        elif method == 'nj':
            if do_bootstrap:
                tree = bootstrap(seqs_list, names, method='nj',
                               n_bootstraps=n_boot, distance_model=dist_model)
            else:
                tree = neighbor_joining(seqs_list, names, distance_model=dist_model)
        elif method == 'mp':
            if do_bootstrap:
                tree = bootstrap(seqs_list, names, method='mp',
                               n_bootstraps=n_boot, distance_model=dist_model)
            else:
                tree = maximum_parsimony(seqs_list, names)
        else:
            tree = upgma(seqs_list, names, distance_model=dist_model)
        
        svg = tree_to_svg(tree, width=600, height=400)
        
        bootstrap_values = []
        def collect_bootstrap(node):
            if not node.is_leaf() and node.bootstrap is not None:
                bootstrap_values.append(round(float(node.bootstrap), 1))
            for child in node.children:
                collect_bootstrap(child)
        if tree.root:
            collect_bootstrap(tree.root)
        
        result['tree'] = {
            'newick': tree.to_newick(),
            'leaves': tree.get_leaves(),
            'svg': svg,
            'method': method,
            'distance_model': dist_model,
            'bootstrap': do_bootstrap,
            'n_bootstrap': n_boot if do_bootstrap else 0,
            'bootstrap_values': bootstrap_values,
        }
    
    elif analysis == 'motif':
        seq = seqs[0]
        motif_type = data.get('motif_type', 'regex')
        pattern = data.get('pattern', '')
        prosite_id = data.get('prosite_id', '')
        
        if motif_type == 'regex' and pattern:
            matches = find_motif_regex(seq, pattern)
            result['motifs'] = matches
        elif motif_type == 'prosite':
            if prosite_id:
                matches = find_prosite_motif(seq, prosite_id=prosite_id)
            else:
                from bio_kit.motif import find_all_prosite_motifs
                matches = find_all_prosite_motifs(seq)
            result['motifs'] = matches[:50]
            result['available_prosite'] = list(PROSITE_PATTERNS.keys())[:20]
    
    elif analysis == 'statistics':
        seq = seqs[0]
        stats = sequence_stats(seq)
        stats['n50'] = n50(seqs)
        
        kmer = kmer_frequency(seq, k=3)
        sorted_kmers = sorted(kmer.items(), key=lambda x: x[1], reverse=True)[:20]
        stats['top_kmers'] = sorted_kmers
        
        stats['rscu'] = rscu_analysis(seq)
        
        result['statistics'] = stats
    
    elif analysis == 'gene_prediction':
        seq = seqs[0]
        promoters = predict_promoters(seq)
        terminators = predict_terminators(seq)
        rbs = predict_rbs(seq)
        
        result['gene_prediction'] = {
            'promoters': promoters[:20],
            'terminators': terminators[:20],
            'rbs': rbs[:20],
        }
    
    elif analysis == 'visualization':
        viz_type = data.get('viz_type', 'gc')
        
        if viz_type == 'gc':
            seq = seqs[0]
            svg = gc_plot(seq, window=50, width=700, height=250)
            result['visualization'] = {'svg': svg, 'type': 'gc'}
        
        elif viz_type == 'alignment':
            if len(seqs) < 2:
                return jsonify({'error': 'Need at least 2 sequences'}), 400
            seqs_list = [s.seq for s in seqs]
            names = [s.id or f'seq{i+1}' for i, s in enumerate(seqs)]
            msa = clustalw_style_msa(seqs_list)
            svg = alignment_colors(msa.aligned_sequences, msa.names, width=700)
            result['visualization'] = {'svg': svg, 'type': 'alignment', 'msa': msa}
        
        elif viz_type == 'heatmap':
            if len(seqs) < 2:
                return jsonify({'error': 'Need at least 2 sequences'}), 400
            seqs_list = [s.seq for s in seqs]
            names = [s.id or f'seq{i+1}' for i, s in enumerate(seqs)]
            svg = similarity_heatmap(seqs_list, names, width=500, height=400)
            result['visualization'] = {'svg': svg, 'type': 'heatmap'}
        
        elif viz_type == 'tree':
            if len(seqs) < 3:
                return jsonify({'error': 'Need at least 3 sequences'}), 400
            seqs_list = [s.seq for s in seqs]
            names = [s.id or f'seq{i+1}' for i, s in enumerate(seqs)]
            tree = upgma(seqs_list, names)
            svg = tree_to_svg(tree, width=600, height=400)
            result['visualization'] = {'svg': svg, 'type': 'tree', 'newick': tree.to_newick()}
    
    return jsonify(result)


@app.route('/api/export', methods=['POST'])
def api_export():
    data = request.json
    content = data.get('content', '')
    fmt = data.get('format', 'fasta')
    
    tmp = tempfile.NamedTemporaryFile(mode='w', suffix=f'.{fmt}', delete=False)
    tmp.write(content)
    tmp.close()
    
    try:
        return send_file(tmp.name, as_attachment=True, download_name=f'bio_kit_result.{fmt}')
    finally:
        pass


if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)
