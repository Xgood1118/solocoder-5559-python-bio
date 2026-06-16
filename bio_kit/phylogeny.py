import random
import math
from copy import deepcopy

class DistanceModel:
    @staticmethod
    def jukes_cantor(p_distance):
        if p_distance >= 0.75:
            return float('inf')
        if p_distance <= 0:
            return 0.0
        return -3/4 * math.log(1 - 4/3 * p_distance)
    
    @staticmethod
    def kimura2(p_distance, transition_ratio=0.5):
        if p_distance <= 0:
            return 0.0
        p = p_distance * transition_ratio
        q = p_distance * (1 - transition_ratio)
        denom = 1 - 2 * p - q
        if denom <= 0:
            return float('inf')
        return -0.5 * math.log(denom) - 0.25 * math.log(1 - 2 * q)
    
    @staticmethod
    def raw_distance(seq1, seq2):
        if len(seq1) != len(seq2):
            min_len = min(len(seq1), len(seq2))
            seq1 = seq1[:min_len]
            seq2 = seq2[:min_len]
        diff = sum(1 for a, b in zip(seq1, seq2) if a != '-' and b != '-' and a != b)
        total = sum(1 for a, b in zip(seq1, seq2) if a != '-' or b != '-')
        return diff / max(total, 1)


class TreeNode:
    def __init__(self, name='', children=None, distance=0.0):
        self.name = name
        self.children = children or []
        self.distance = distance
        self.bootstrap = None
    
    def is_leaf(self):
        return len(self.children) == 0
    
    def get_leaves(self):
        if self.is_leaf():
            return [self.name]
        leaves = []
        for child in self.children:
            leaves.extend(child.get_leaves())
        return leaves
    
    def to_newick(self):
        if self.is_leaf():
            return f"{self.name}:{self.distance:.6f}"
        child_strs = []
        for child in self.children:
            child_strs.append(child.to_newick())
        inner = ','.join(child_strs)
        if self.bootstrap is not None:
            return f"({inner}){self.bootstrap:.0f}:{self.distance:.6f}"
        return f"({inner}):{self.distance:.6f}"


class Tree:
    def __init__(self, root=None):
        self.root = root
    
    def to_newick(self):
        if not self.root:
            return ';'
        return self.root.to_newick() + ';'
    
    def get_leaves(self):
        if not self.root:
            return []
        return self.root.get_leaves()
    
    def __repr__(self):
        return f"Tree(leaves={len(self.get_leaves())})"


def _p_distance_matrix(sequences):
    n = len(sequences)
    dist = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            d = DistanceModel.raw_distance(sequences[i], sequences[j])
            dist[i][j] = d
            dist[j][i] = d
    return dist


def _apply_distance_model(dist_matrix, model='jukes_cantor'):
    n = len(dist_matrix)
    result = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j:
                result[i][j] = 0.0
            else:
                if model == 'jukes_cantor':
                    result[i][j] = DistanceModel.jukes_cantor(dist_matrix[i][j])
                elif model == 'kimura2':
                    result[i][j] = DistanceModel.kimura2(dist_matrix[i][j])
                else:
                    result[i][j] = dist_matrix[i][j]
    return result


def upgma(sequences, names=None, distance_model='jukes_cantor'):
    if isinstance(sequences[0], str):
        seqs = [s.upper() for s in sequences]
    else:
        seqs = [s.seq.upper() for s in sequences]
        if names is None:
            names = [s.id for s in sequences]
    
    if names is None:
        names = [f'seq{i+1}' for i in range(len(seqs))]
    
    n = len(seqs)
    p_dist = _p_distance_matrix(seqs)
    dist = _apply_distance_model(p_dist, distance_model)
    
    clusters = [[i] for i in range(n)]
    nodes = [TreeNode(name=names[i]) for i in range(n)]
    heights = [0.0] * n
    
    while len(clusters) > 1:
        min_dist = float('inf')
        best_i, best_j = 0, 1
        
        for i in range(len(clusters)):
            for j in range(i + 1, len(clusters)):
                total = 0.0
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
        new_height = min_dist / 2.0
        
        left_node = nodes[best_i]
        right_node = nodes[best_j]
        left_node.distance = new_height - heights[best_i]
        right_node.distance = new_height - heights[best_j]
        
        new_node = TreeNode(children=[left_node, right_node])
        
        if best_i < best_j:
            clusters.pop(best_j)
            clusters.pop(best_i)
            nodes.pop(best_j)
            nodes.pop(best_i)
            heights.pop(best_j)
            heights.pop(best_i)
        else:
            clusters.pop(best_i)
            clusters.pop(best_j)
            nodes.pop(best_i)
            nodes.pop(best_j)
            heights.pop(best_i)
            heights.pop(best_j)
        
        clusters.insert(best_i, new_cluster)
        nodes.insert(best_i, new_node)
        heights.insert(best_i, new_height)
    
    return Tree(root=nodes[0])


def neighbor_joining(sequences, names=None, distance_model='jukes_cantor'):
    if isinstance(sequences[0], str):
        seqs = [s.upper() for s in sequences]
    else:
        seqs = [s.seq.upper() for s in sequences]
        if names is None:
            names = [s.id for s in sequences]
    
    if names is None:
        names = [f'seq{i+1}' for i in range(len(seqs))]
    
    n = len(seqs)
    p_dist = _p_distance_matrix(seqs)
    dist = _apply_distance_model(p_dist, distance_model)
    
    active = list(range(n))
    nodes = [TreeNode(name=names[i]) for i in range(n)]
    
    while len(active) > 2:
        m = len(active)
        
        r = [0.0] * m
        for i in range(m):
            r[i] = sum(dist[active[i]][active[j]] for j in range(m)) / (m - 2)
        
        min_q = float('inf')
        best_i, best_j = 0, 1
        
        for i in range(m):
            for j in range(i + 1, m):
                q = dist[active[i]][active[j]] - r[i] - r[j]
                if q < min_q:
                    min_q = q
                    best_i, best_j = i, j
        
        idx_i = active[best_i]
        idx_j = active[best_j]
        
        d_i = (dist[idx_i][idx_j] + r[best_i] - r[best_j]) / 2
        d_j = (dist[idx_i][idx_j] + r[best_j] - r[best_i]) / 2
        
        left_node = nodes[idx_i]
        right_node = nodes[idx_j]
        left_node.distance = max(d_i, 0.001)
        right_node.distance = max(d_j, 0.001)
        
        new_node = TreeNode(children=[left_node, right_node])
        new_idx = max(active) + 1
        
        new_dists = {}
        for k_idx in active:
            if k_idx == idx_i or k_idx == idx_j:
                continue
            d = (dist[idx_i][k_idx] + dist[idx_j][k_idx] - dist[idx_i][idx_j]) / 2
            new_dists[k_idx] = max(d, 0.0)
        
        for row in dist:
            row.append(0.0)
        dist.append([0.0] * (len(dist[0])))
        
        for k_idx in active:
            if k_idx == idx_i or k_idx == idx_j:
                continue
            dist[new_idx][k_idx] = new_dists[k_idx]
            dist[k_idx][new_idx] = new_dists[k_idx]
        
        nodes.append(new_node)
        
        remove_i = best_i if best_i > best_j else best_j
        remove_j = best_j if best_i > best_j else best_i
        active.pop(remove_i)
        active.pop(remove_j)
        active.append(new_idx)
    
    if len(active) == 2:
        idx_i, idx_j = active[0], active[1]
        d = dist[idx_i][idx_j] / 2
        nodes[idx_i].distance = max(d, 0.001)
        nodes[idx_j].distance = max(d, 0.001)
        root = TreeNode(children=[nodes[idx_i], nodes[idx_j]])
    else:
        root = nodes[active[0]]
    
    return Tree(root=root)


def maximum_parsimony(sequences, names=None):
    if isinstance(sequences[0], str):
        seqs_list = [s.upper() for s in sequences]
    else:
        seqs_list = [s.seq.upper() for s in sequences]
        if names is None:
            names = [s.id for s in sequences]
    
    if names is None:
        names = [f'seq{i+1}' for i in range(len(seqs_list))]
    
    n = len(seqs_list)
    if n < 3:
        return upgma(seqs_list, names)
    
    seq_len = min(len(s) for s in seqs_list)
    seq_chars = [list(s[:seq_len]) for s in seqs_list]
    
    def fitch_parsimony(tree):
        total_score = 0
        
        def postorder(node):
            nonlocal total_score
            if node.is_leaf():
                idx = names.index(node.name) if node.name in names else 0
                if idx >= len(seq_chars):
                    idx = 0
                return set()
            
            child_sets = []
            for child in node.children:
                child_sets.append(postorder(child))
            
            if not child_sets:
                return set()
            
            intersection = set(child_sets[0])
            for s in child_sets[1:]:
                intersection &= s
            
            if intersection:
                return intersection
            else:
                total_score += 1
                union = set()
                for s in child_sets:
                    union |= s
                return union
        
        for site in range(seq_len):
            def postorder_site(node):
                nonlocal total_score
                if node.is_leaf():
                    idx = names.index(node.name) if node.name in names else 0
                    if idx >= len(seq_chars):
                        idx = 0
                    return {seq_chars[idx][site]}
                
                child_sets = []
                for child in node.children:
                    child_sets.append(postorder_site(child))
                
                if not child_sets:
                    return set()
                
                intersection = set(child_sets[0])
                for s in child_sets[1:]:
                    intersection &= s
                
                if intersection:
                    return intersection
                else:
                    total_score += 1
                    union = set()
                    for s in child_sets:
                        union |= s
                    return union
            
            postorder_site(tree.root)
        
        return total_score
    
    tree = upgma(seqs_list, names)
    
    best_tree = deepcopy(tree)
    best_score = fitch_parsimony(best_tree)
    
    def get_internal_nodes(node):
        nodes = []
        if not node.is_leaf():
            nodes.append(node)
            for child in node.children:
                nodes.extend(get_internal_nodes(child))
        return nodes
    
    improved = True
    iterations = 0
    while improved and iterations < 5:
        improved = False
        iterations += 1
        
        internal_nodes = get_internal_nodes(best_tree.root)
        
        for node in internal_nodes:
            if len(node.children) < 2:
                continue
            
            children = node.children
            
            for i in range(len(children)):
                for j in range(i + 1, len(children)):
                    children[i], children[j] = children[j], children[i]
                    
                    score = fitch_parsimony(best_tree)
                    if score < best_score:
                        best_score = score
                        best_tree = deepcopy(best_tree)
                        improved = True
                    else:
                        children[i], children[j] = children[j], children[i]
    
    return best_tree


def bootstrap(sequences, names=None, method='upgma', n_bootstraps=100, distance_model='jukes_cantor'):
    if isinstance(sequences[0], str):
        seqs = [s.upper() for s in sequences]
    else:
        seqs = [s.seq.upper() for s in sequences]
        if names is None:
            names = [s.id for s in sequences]
    
    if names is None:
        names = [f'seq{i+1}' for i in range(len(seqs))]
    
    if method == 'upgma':
        tree_func = lambda s, n: upgma(s, n, distance_model)
    elif method == 'nj':
        tree_func = lambda s, n: neighbor_joining(s, n, distance_model)
    elif method == 'mp':
        tree_func = lambda s, n: maximum_parsimony(s, n)
    else:
        tree_func = lambda s, n: upgma(s, n, distance_model)
    
    ref_tree = tree_func(seqs, names)
    
    seq_len = min(len(s) for s in seqs)
    
    bipartition_counts = {}
    
    def get_bipartitions(tree):
        bipartitions = []
        
        def traverse(node):
            if node.is_leaf():
                return {node.name}
            
            left_set = traverse(node.children[0])
            right_set = traverse(node.children[1])
            
            all_leaves = left_set | right_set
            left_sorted = tuple(sorted(left_set))
            right_sorted = tuple(sorted(right_set))
            
            if left_sorted < right_sorted:
                bipartitions.append((left_sorted, right_sorted))
            else:
                bipartitions.append((right_sorted, left_sorted))
            
            return all_leaves
        
        traverse(tree.root)
        return set(bipartitions)
    
    ref_bipartitions = get_bipartitions(ref_tree)
    for bp in ref_bipartitions:
        bipartition_counts[bp] = 0
    
    for _ in range(n_bootstraps):
        indices = [random.randint(0, seq_len - 1) for _ in range(seq_len)]
        boot_seqs = []
        for seq in seqs:
            boot_seq = ''.join(seq[i] for i in indices)
            boot_seqs.append(boot_seq)
        
        try:
            boot_tree = tree_func(boot_seqs, names)
            boot_bipartitions = get_bipartitions(boot_tree)
            for bp in ref_bipartitions:
                if bp in boot_bipartitions:
                    bipartition_counts[bp] += 1
        except:
            continue
    
    def assign_bootstrap(node):
        if node.is_leaf():
            return {node.name}
        
        left_set = assign_bootstrap(node.children[0])
        right_set = assign_bootstrap(node.children[1])
        
        left_sorted = tuple(sorted(left_set))
        right_sorted = tuple(sorted(right_set))
        
        if left_sorted < right_sorted:
            bp = (left_sorted, right_sorted)
        else:
            bp = (right_sorted, left_sorted)
        
        if bp in bipartition_counts:
            node.bootstrap = (bipartition_counts[bp] / n_bootstraps) * 100
        
        return left_set | right_set
    
    assign_bootstrap(ref_tree.root)
    
    return ref_tree
