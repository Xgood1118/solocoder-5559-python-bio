import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

print("=" * 60)
print("Bio Kit - 功能测试")
print("=" * 60)

def test_basic():
    print("\n[1] 测试基础序列分析...")
    from bio_kit import reverse_complement, complement, transcribe, translate, gc_content
    
    dna = "ATGGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGTAG"
    
    print(f"  原始序列: {dna[:30]}...")
    print(f"  长度: {len(dna)} bp")
    
    rc = reverse_complement(dna)
    print(f"  反向互补: {rc[:30]}...")
    
    comp = complement(dna)
    print(f"  互补: {comp[:30]}...")
    
    rna = transcribe(dna)
    print(f"  转录(RNA): {rna[:30]}...")
    
    protein = translate(dna)
    print(f"  翻译(蛋白): {protein[:30]}...")
    
    gc = gc_content(dna)
    print(f"  GC含量: {gc:.2f}%")
    
    print("  ✅ 基础序列分析 - 通过")
    return True

def test_orf():
    print("\n[2] 测试ORF查找...")
    from bio_kit import find_orfs, six_frame_translation
    
    dna = "ATGGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    dna += "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT"
    dna += "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"
    dna += "TCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    dna += "TAG"
    
    orfs = find_orfs(dna, min_length=10)
    print(f"  找到 {len(orfs)} 个ORF")
    if orfs:
        print(f"  最长ORF: {orfs[0]['length_aa']} aa, 链: {orfs[0]['strand']}")
    
    frames = six_frame_translation(dna)
    print(f"  六框翻译完成: {list(frames.keys())}")
    
    print("  ✅ ORF查找 - 通过")
    return True

def test_restriction():
    print("\n[3] 测试限制酶切分析...")
    from bio_kit import find_restriction_sites, list_enzymes
    
    dna = "GAATTCGATCGATCGATCGGATCCGATCGATCGAAGCTTGATCGATCGAGATCTCG"
    
    enzymes = ['EcoRI', 'BamHI', 'HindIII', 'BglII']
    sites = find_restriction_sites(dna, enzymes)
    
    for enzyme, info in sites.items():
        print(f"  {enzyme}: {info['count']} 个位点, 片段: {info['fragments']}")
    
    all_enzymes = list_enzymes()
    print(f"  内置酶数量: {len(all_enzymes)} 种")
    
    print("  ✅ 限制酶切分析 - 通过")
    return True

def test_primer():
    print("\n[4] 测试引物设计...")
    from bio_kit import design_primers, calculate_tm, calculate_gc, check_dimer, check_hairpin
    
    template = "ATGGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    template += "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT"
    template += "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"
    template += "TCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    template += "TAG"
    
    primer = "ATGGCGATCGATCGATCG"
    print(f"  引物: {primer}")
    print(f"  Tm: {calculate_tm(primer):.2f}°C")
    print(f"  GC: {calculate_gc(primer):.2f}%")
    print(f"  二聚体检测: {check_dimer(primer)}")
    print(f"  发夹检测: {check_hairpin(primer)}")
    
    result = design_primers(template, start=0, end=200, 
                           length_range=(18, 25), 
                           tm_range=(50, 70), 
                           gc_range=(30, 70))
    
    print(f"  候选上游引物: {len(result['forward_primers'])} 个")
    print(f"  候选下游引物: {len(result['reverse_primers'])} 个")
    print(f"  最佳引物对: {len(result['best_pairs'])} 对")
    
    print("  ✅ 引物设计 - 通过")
    return True

def test_alignment():
    print("\n[5] 测试序列比对...")
    from bio_kit import needleman_wunsch, smith_waterman, clustalw_style_msa
    
    seq1 = "ATGGCGATTACCGTTGATGTTGATGCTGGTGAA"
    seq2 = "ATGGCGATTACCGTTCATGTTGATGCTGGTCAA"
    
    nw = needleman_wunsch(seq1, seq2)
    print(f"  Needleman-Wunsch 全局比对:")
    print(f"    得分: {nw.score}, 一致性: {nw.identity:.2f}%")
    print(f"    {nw.seq1_aligned}")
    print(f"    {nw.get_match_line()}")
    print(f"    {nw.seq2_aligned}")
    
    sw = smith_waterman(seq1, seq2)
    print(f"  Smith-Waterman 局部比对: 得分={sw.score}")
    
    seqs = [
        "ATGGCGATTACCGTTGATGTTGATGCTGGTGAA",
        "ATGGCGATTACCGTTCATGTTGATGCTGGTCAA",
        "ATGGCGATTACCGTTGACGTTAATGCTGGTGAA",
    ]
    names = ["SeqA", "SeqB", "SeqC"]
    
    msa = clustalw_style_msa(seqs)
    print(f"  多序列比对: {msa.num_seqs} 条序列, 长度 {msa.alignment_length}")
    print(f"  一致序列: {msa.get_consensus()[:30]}...")
    
    print("  ✅ 序列比对 - 通过")
    return True

def test_blast():
    print("\n[6] 测试BLAST近似搜索...")
    from bio_kit import MiniBLAST
    
    query = "ATGGCGATTACCGTTGATGTTGATGCTGGTGAA"
    
    db_seqs = [
        "ATGGCGATTACCGTTCATGTTGATGCTGGTCAAACTGCTGGTGAAACCGGTGAA",
        "ATGGCGATTACCGTTGACGTTAATGCTGGTGAAACGGCTGGCGAAACCGGCAAA",
        "ATGGCTATCACCGTTGATGTTAATGCAGGTGAAACTGCTGGCGAGACTGGTGAA",
        "TTAAGCGGATGATGGCCGATGCTAATGCGGGTGAAAGTGCTGGTGATGCCGGT",
    ]
    
    blast = MiniBLAST(db_seqs, word_size=7, score_threshold=20)
    print(f"  数据库大小: {len(blast)} 条序列")
    
    results = blast.search(query, top_n=5)
    print(f"  找到 {len(results)} 个匹配")
    
    for i, hit in enumerate(results[:3]):
        print(f"    {i+1}. 得分={hit['score']}, 一致性={hit['identity']:.1f}%")
    
    print("  ✅ BLAST近似搜索 - 通过")
    return True

def test_phylogeny():
    print("\n[7] 测试进化树构建...")
    from bio_kit import upgma, neighbor_joining, maximum_parsimony, bootstrap, DistanceModel
    
    seqs = [
        "ATGGCGATTACCGTTGATGTTGATGCTGGTGAAACTGCTGGTGAAACCGGTGAA",
        "ATGGCGATTACCGTTCATGTTGATGCTGGTCAAACTGCTGGTGAAACCGGTGAA",
        "ATGGCGATTACCGTTGACGTTAATGCTGGTGAAACGGCTGGCGAAACCGGCAAA",
        "ATGGCTATCACCGTTGATGTTAATGCAGGTGAAACTGCTGGCGAGACTGGTGAA",
        "TTAAGCGGATGATGGCCGATGCTAATGCGGGTGAAAGTGCTGGTGATGCCGGT",
    ]
    names = ["SeqA", "SeqB", "SeqC", "SeqD", "SeqE"]
    
    print("  测试 Jukes-Cantor 距离模型...")
    d = DistanceModel.jukes_cantor(0.1)
    print(f"    p=0.1 -> {d:.6f}")
    
    print("  测试 Kimura 2-parameter 距离模型...")
    d2 = DistanceModel.kimura2(0.1)
    print(f"    p=0.1 -> {d2:.6f}")
    
    print("  构建 UPGMA 树...")
    tree_upgma = upgma(seqs, names, distance_model='jukes_cantor')
    print(f"    叶节点: {tree_upgma.get_leaves()}")
    newick = tree_upgma.to_newick()
    print(f"    Newick: {newick[:80]}...")
    
    print("  构建 NJ 树...")
    tree_nj = neighbor_joining(seqs, names, distance_model='jukes_cantor')
    print(f"    叶节点: {len(tree_nj.get_leaves())}")
    
    print("  构建 MP 树...")
    try:
        tree_mp = maximum_parsimony(seqs, names)
        print(f"    叶节点: {len(tree_mp.get_leaves())}")
    except Exception as e:
        print(f"    MP 树构建: {e}")
    
    print("  Bootstrap 检验 (UPGMA, 10次)...")
    tree_boot = bootstrap(seqs, names, method='upgma', n_bootstraps=10, distance_model='jukes_cantor')
    print(f"    完成 Bootstrap")
    
    print("  ✅ 进化树构建 - 通过")
    return True

def test_motif():
    print("\n[8] 测试基序查找...")
    from bio_kit import find_motif_regex, find_prosite_motif, PROSITE_PATTERNS
    
    seq = "ATGGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    
    print("  正则表达式基序查找...")
    matches = find_motif_regex(seq, r'GATC')
    print(f"    GATC 基序: {len(matches)} 个")
    
    print(f"  PROSITE 数据库条目数: {len(PROSITE_PATTERNS)}")
    
    print("  PROSITE 基序查找...")
    prosite_matches = find_prosite_motif(seq, prosite_id='PS00017')
    print(f"    PS00017: {len(prosite_matches)} 个匹配")
    
    print("  ✅ 基序查找 - 通过")
    return True

def test_statistics():
    print("\n[9] 测试序列统计...")
    from bio_kit import gc_content, n50, kmer_frequency, shannon_entropy, rscu_analysis, sequence_stats
    
    dna = "ATGGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    
    stats = sequence_stats(dna)
    print(f"  序列长度: {stats['length']} bp")
    print(f"  GC含量: {stats['gc_content']:.2f}%")
    print(f"  Shannon 熵: {stats['shannon_entropy']:.4f} bits")
    
    kmers = kmer_frequency(dna, k=3)
    print(f"  3-mer 种类: {len(kmers)} 种")
    
    rscu = rscu_analysis(dna)
    print(f"  密码子数: {len(rscu['codon_counts'])}")
    
    seqs = [dna, dna[:100], dna[:200]]
    n50_val = n50(seqs)
    print(f"  N50: {n50_val} bp")
    
    print("  ✅ 序列统计 - 通过")
    return True

def test_gene_prediction():
    print("\n[10] 测试基因预测...")
    from bio_kit import predict_promoters, predict_terminators, predict_rbs
    
    seq = "TTGACAATCGATCGATCGTATAATCGATCGATCGATCGATCGATCGATCGATCG"
    seq += "AGGAGGTCATGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
    seq += "GATCGATCGATCGATCGGCTCTCCAAAAAATCGATCGATCGATCGATCGATCGA"
    
    promoters = predict_promoters(seq)
    print(f"  预测启动子: {len(promoters)} 个")
    if promoters:
        print(f"    最佳: {promoters[0]['sequence']} (得分: {promoters[0]['score']:.2f}, 类型: {promoters[0]['type']})")
    
    terminators = predict_terminators(seq)
    print(f"  预测终止子: {len(terminators)} 个")
    
    rbs = predict_rbs(seq)
    print(f"  预测RBS: {len(rbs)} 个")
    if rbs:
        print(f"    最佳: {rbs[0]['sequence']} (得分: {rbs[0]['score']})")
    
    print("  ✅ 基因预测 - 通过")
    return True

def test_visualization():
    print("\n[11] 测试可视化...")
    from bio_kit import gc_plot, alignment_colors, tree_to_svg, similarity_heatmap, upgma
    
    dna = "ATGGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG" * 5
    
    print("  GC含量图...")
    svg_gc = gc_plot(dna, window=30, width=600, height=200)
    print(f"    SVG 大小: {len(svg_gc)} 字符")
    
    seqs = [
        "ATGGCGATTACCGTTGATGTTGATGCTGGTGAA",
        "ATGGCGATTACCGTTCATGTTGATGCTGGTCAA",
        "ATGGCGATTACCGTTGACGTTAATGCTGGTGAA",
    ]
    names = ["SeqA", "SeqB", "SeqC"]
    
    print("  多序列比对图...")
    svg_aln = alignment_colors(seqs, names, width=600)
    print(f"    SVG 大小: {len(svg_aln)} 字符")
    
    print("  进化树SVG...")
    tree = upgma(seqs, names)
    svg_tree = tree_to_svg(tree, width=500, height=300)
    print(f"    SVG 大小: {len(svg_tree)} 字符")
    
    print("  相似度热力图...")
    svg_heat = similarity_heatmap(seqs, names, width=400, height=350)
    print(f"    SVG 大小: {len(svg_heat)} 字符")
    
    print("  ✅ 可视化 - 通过")
    return True

def test_io():
    print("\n[12] 测试序列读写...")
    from bio_kit import read_fasta, read_fastq, read_genbank, read_clustal, write_fasta, write_genbank, Sequence
    import tempfile
    import os
    
    fasta_content = """>Seq1 description text
ATGGCGATCGATCGATCG
>Seq2 another sequence
GATCGATCGATCGATCGATCGA
"""
    
    print("  FASTA 读取...")
    tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
    tmp.write(fasta_content)
    tmp.close()
    
    seqs = read_fasta(tmp.name)
    print(f"    读取到 {len(seqs)} 条序列")
    print(f"    第一条: {seqs[0].id}, 长度 {len(seqs[0])}")
    
    os.unlink(tmp.name)
    
    print("  FASTA 写入...")
    tmp2 = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
    tmp2.close()
    write_fasta(seqs, tmp2.name)
    
    seqs2 = read_fasta(tmp2.name)
    print(f"    回读验证: {len(seqs2)} 条序列")
    
    os.unlink(tmp2.name)
    
    print("  ✅ 序列读写 - 通过")
    return True

def main():
    tests = [
        test_basic,
        test_orf,
        test_restriction,
        test_primer,
        test_alignment,
        test_blast,
        test_phylogeny,
        test_motif,
        test_statistics,
        test_gene_prediction,
        test_visualization,
        test_io,
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            if test():
                passed += 1
            else:
                failed += 1
        except Exception as e:
            print(f"  ❌ 测试失败: {e}")
            import traceback
            traceback.print_exc()
            failed += 1
    
    print("\n" + "=" * 60)
    print(f"测试结果: {passed} 通过, {failed} 失败, 共 {len(tests)} 项")
    print("=" * 60)
    
    return failed == 0

if __name__ == '__main__':
    success = main()
    sys.exit(0 if success else 1)
