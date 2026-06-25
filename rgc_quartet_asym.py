#!/usr/bin/env python3
"""
rgc_quartet_asym.py - Quartet Asymmetry Test for Binary Characters

Implements the quartet asymmetry test from Springer et al. (2020)
DOI: 10.1093/jhered/esz076

Given a binary character matrix (RGC/indel data) and a species tree, this program:
1. Enumerates all quadripartitions defined by internal branches of the tree.
2. For each quadripartition, enumerates only VALID quartets: one taxon drawn from
   each of the four sub-groups (a_taxa, b_taxa, c_taxa, d_taxa).  This guarantees
   that the tree-concordant topology is always ab|cd by construction and prevents
   quartets with two taxa from the same sub-group being mis-attributed to a
   higher-level branch.
3. Tests whether the two minority quartets have symmetric counts (binomial test).
4. Applies multiple test corrections (BH, BY, Bonferroni, Holm-Bonferroni)
   within each quadripartition independently.

The majority quartet (ab|cd) should match the species tree topology.
Under ILS without introgression, the two minority quartets (ac|bd and ad|bc)
should have equal probabilities, so asymmetry suggests gene flow.
"""

import argparse
import sys
import re
import os
from datetime import datetime
from math import comb, log
from collections import defaultdict


def timestamp():
    """Return current timestamp string."""
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def log_progress(message):
    """Print timestamped progress message."""
    print("[{}] {}".format(timestamp(), message))


def parse_args():
    parser = argparse.ArgumentParser(
        description='Quartet Asymmetry Test for Binary Characters (RGC/Indel data)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python rgc_quartet_asym.py -i data.phy -t species.tre -o results
  python rgc_quartet_asym.py -i data.nex -t tree.nwk -o results --alpha 0.01
  python rgc_quartet_asym.py -i data.phy -t tree.nwk -o results --correction BH BY

Output files:
  <output>.quartet_asym.tsv      - Detailed per-quartet results
  <output>.quadripartition_summary.tsv - Summary statistics per quadripartition
  <output>.labeled_tree.tre      - Species tree with internal node labels
        """
    )
    parser.add_argument('--input', '-i', required=True,
                        help='Input character matrix (PHYLIP or NEXUS format)')
    parser.add_argument('--tree', '-t', required=True,
                        help='Species tree in Newick format')
    parser.add_argument('--output', '-o', required=True,
                        help='Base name for output files')
    parser.add_argument('--alpha', type=float, default=0.05,
                        help='Significance level (default: 0.05)')
    parser.add_argument('--correction', nargs='+', 
                        choices=['BH', 'BY', 'bonferroni', 'holm'],
                        default=['BH', 'BY', 'bonferroni', 'holm'],
                        help='Multiple test correction methods (default: BH BY bonferroni holm)')
    return parser.parse_args()


def detect_format(filename):
    """Detect whether input file is NEXUS or PHYLIP format."""
    with open(filename, 'r') as f:
        first_line = f.readline().strip()
        if first_line.upper().startswith('#NEXUS'):
            return 'nexus'
        else:
            return 'phylip'


def parse_phylip(filename):
    """Parse relaxed PHYLIP format (taxon names of any length)."""
    taxa = []
    sequences = {}
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # First line: ntax nchar
    header = lines[0].strip().split()
    ntax = int(header[0])
    nchar = int(header[1])
    
    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue
        
        # Split on first whitespace block
        parts = line.split(None, 1)
        if len(parts) != 2:
            continue
        
        taxon = parts[0]
        seq = parts[1].replace(' ', '').replace('\t', '')
        
        taxa.append(taxon)
        sequences[taxon] = seq
    
    if len(taxa) != ntax:
        log_progress("Warning: Expected {} taxa, found {}".format(ntax, len(taxa)))
    
    if taxa and len(sequences[taxa[0]]) != nchar:
        log_progress("Warning: Expected {} characters, found {}".format(nchar, len(sequences[taxa[0]])))
    
    return taxa, sequences


def parse_nexus(filename):
    """Parse NEXUS format character matrix."""
    taxa = []
    sequences = {}
    
    with open(filename, 'r') as f:
        content = f.read()
    
    # Find matrix block
    matrix_match = re.search(r'matrix\s+(.*?)\s*;\s*end;', content, 
                             re.IGNORECASE | re.DOTALL)
    if not matrix_match:
        raise ValueError("Could not find MATRIX block in NEXUS file")
    
    matrix_text = matrix_match.group(1)
    
    for line in matrix_text.strip().split('\n'):
        line = line.strip()
        if not line or line.startswith('[') or line.startswith('#'):
            continue
        
        parts = line.split(None, 1)
        if len(parts) != 2:
            continue
        
        taxon = parts[0]
        seq = parts[1].replace(' ', '').replace('\t', '')
        
        # Remove any bracketed comments within sequence
        seq = re.sub(r'\[[^\]]*\]', '', seq)
        
        taxa.append(taxon)
        sequences[taxon] = seq
    
    return taxa, sequences


# =============================================================================
# Newick Tree Parsing and Manipulation
# =============================================================================

class TreeNode:
    """Simple tree node class for phylogenetic trees."""
    
    def __init__(self, name=None):
        self.name = name
        self.children = []
        self.parent = None
        self.branch_length = None
        self.node_id = None  # Numeric ID for internal nodes
    
    def is_leaf(self):
        return len(self.children) == 0
    
    def get_leaves(self):
        """Return set of leaf names in subtree rooted at this node."""
        if self.is_leaf():
            return {self.name}
        leaves = set()
        for child in self.children:
            leaves.update(child.get_leaves())
        return leaves
    
    def get_all_nodes(self):
        """Return list of all nodes in subtree (preorder)."""
        nodes = [self]
        for child in self.children:
            nodes.extend(child.get_all_nodes())
        return nodes


def parse_newick(newick_str):
    """
    Parse a Newick string into a tree structure.
    Returns the root node.
    """
    newick_str = newick_str.strip()
    if newick_str.endswith(';'):
        newick_str = newick_str[:-1]
    
    def parse_node(s, pos):
        """Parse a node starting at position pos, return (node, new_pos)."""
        node = TreeNode()
        
        if pos < len(s) and s[pos] == '(':
            # Internal node with children
            pos += 1  # Skip '('
            
            while True:
                child, pos = parse_node(s, pos)
                child.parent = node
                node.children.append(child)
                
                if pos >= len(s):
                    break
                if s[pos] == ',':
                    pos += 1  # Skip ','
                elif s[pos] == ')':
                    pos += 1  # Skip ')'
                    break
                else:
                    raise ValueError(f"Unexpected character at position {pos}: {s[pos]}")
        
        # Parse name (everything up to : or , or ) or ;)
        name_start = pos
        while pos < len(s) and s[pos] not in ':,);':
            pos += 1
        
        if pos > name_start:
            node.name = s[name_start:pos].strip()
        
        # Parse branch length if present
        if pos < len(s) and s[pos] == ':':
            pos += 1  # Skip ':'
            bl_start = pos
            while pos < len(s) and s[pos] not in ',);':
                pos += 1
            if pos > bl_start:
                try:
                    node.branch_length = float(s[bl_start:pos])
                except ValueError:
                    pass
        
        return node, pos
    
    root, _ = parse_node(newick_str, 0)
    return root


def assign_node_ids(root):
    """
    Assign numeric IDs to internal nodes.
    Returns dict mapping node_id -> node.
    """
    node_map = {}
    node_id = 1
    
    for node in root.get_all_nodes():
        if not node.is_leaf():
            node.node_id = node_id
            node_map[node_id] = node
            node_id += 1
    
    return node_map


def tree_to_newick(node, include_node_ids=True, include_branch_lengths=True):
    """Convert tree back to Newick string with optional node IDs."""
    if node.is_leaf():
        result = node.name if node.name else ""
    else:
        child_strs = [tree_to_newick(child, include_node_ids, include_branch_lengths) 
                      for child in node.children]
        result = "(" + ",".join(child_strs) + ")"
        if include_node_ids and node.node_id is not None:
            result += str(node.node_id)
        elif node.name:
            result += node.name
    
    if include_branch_lengths and node.branch_length is not None:
        result += ":" + str(node.branch_length)
    
    return result


# =============================================================================
# Quadripartition Extraction
# =============================================================================

def _enumerate_group_pairs(groups):
    """
    Given a list of groups (sets of taxa), enumerate ways to partition them into two parts.
    Each part can contain one or more original groups combined.
    
    For a binary tree branch, each side has exactly 2 subtrees (or one subtree and 
    "above"), so we just return the groups as-is if there are exactly 2.
    If there are more, we need to combine them in all possible ways.
    """
    if len(groups) == 0:
        return []
    if len(groups) == 1:
        # Single group: can't split further, but we need pairs
        return []
    if len(groups) == 2:
        return [(groups[0], groups[1])]
    
    # More than 2 groups: enumerate all ways to split into two non-empty parts
    # This shouldn't happen for a bifurcating tree, but handle gracefully
    result = []
    n = len(groups)
    for i in range(1, 2**(n-1)):  # All non-trivial bipartitions
        part1 = set()
        part2 = set()
        for j in range(n):
            if i & (1 << j):
                part1.update(groups[j])
            else:
                part2.update(groups[j])
        if part1 and part2:
            result.append((part1, part2))
    
    # Remove duplicates (since bipartition {A,B}|{C,D} = {C,D}|{A,B})
    unique = []
    seen = set()
    for p1, p2 in result:
        key = (frozenset(p1), frozenset(p2))
        key_rev = (frozenset(p2), frozenset(p1))
        if key not in seen and key_rev not in seen:
            unique.append((p1, p2))
            seen.add(key)
    
    return unique


def _format_quadripartition(a_taxa, b_taxa, c_taxa, d_taxa, taxon_order, node1, node2):
    """
    Format a quadripartition with proper ordering.
    
    Ordering rules:
    1. Each group is defined by its smallest taxon
    2. Within each side of the pipe, order by smallest taxon (leftmost has smaller)
    3. The leftmost side overall has the globally smallest taxon
    """
    def min_taxon_num(taxa_set):
        return min(taxon_order[t] for t in taxa_set)
    
    def sort_taxa(taxa_set):
        return sorted(taxa_set, key=lambda t: taxon_order[t])
    
    # Order groups within each side
    side1 = [a_taxa, b_taxa]
    side2 = [c_taxa, d_taxa]
    
    # Sort groups within each side by their minimum taxon
    side1.sort(key=min_taxon_num)
    side2.sort(key=min_taxon_num)
    
    # Determine which side should be on the left (contains the globally smallest taxon)
    min1 = min(min_taxon_num(g) for g in side1)
    min2 = min(min_taxon_num(g) for g in side2)
    
    if min2 < min1:
        side1, side2 = side2, side1
    
    # Now side1 is the left side, side2 is the right side
    # Format each group
    def format_group(taxa_set):
        sorted_taxa = sort_taxa(taxa_set)
        return "(" + ",".join(sorted_taxa) + ")"
    
    left_str = ",".join(format_group(g) for g in side1)
    right_str = ",".join(format_group(g) for g in side2)
    
    formatted = left_str + "|" + right_str
    
    # Assign a, b, c, d labels based on position
    # a = first group left, b = second group left
    # c = first group right, d = second group right
    
    return {
        'node1': node1.node_id,
        'node2': node2.node_id if node2.node_id else 'root',
        'groups': [side1[0], side1[1], side2[0], side2[1]],
        'formatted': formatted,
        'a_taxa': side1[0],
        'b_taxa': side1[1],
        'c_taxa': side2[0],
        'd_taxa': side2[1]
    }


def get_all_quadripartitions(root, taxon_order):
    """
    Get all unique quadripartitions from the tree.
    Each internal branch of a bifurcating tree with n taxa produces one quadripartition.

    For each branch, we identify:
    - The four sub-groups (a_taxa, b_taxa on child side; c_taxa, d_taxa on parent side)
      used to enumerate valid quartets with one taxon per sub-group.
    - child_side / parent_side (complete sets) are retained for reference.
    """
    quadripartitions = []
    all_taxa = root.get_leaves()
    
    # Process each internal node and its relationship to children
    def process_branch(parent_node, child_node):
        """Process the branch between parent_node and child_node."""
        # Taxa on child's side (complete set)
        child_side = child_node.get_leaves()
        # Taxa on parent's side (complete set)
        parent_side = all_taxa - child_side
        
        if len(child_side) < 2 or len(parent_side) < 2:
            return None
        
        # Get the sub-groups on child's side (for display purposes)
        if child_node.is_leaf():
            child_groups = [{child_node.name}]
        else:
            child_groups = [c.get_leaves() for c in child_node.children]
        
        # Get the sub-groups on parent's side (for display purposes)
        parent_groups = []
        for sibling in parent_node.children:
            if sibling != child_node:
                parent_groups.append(sibling.get_leaves())
        # Add taxa "above" parent if it has a parent
        if parent_node.parent:
            above = all_taxa - parent_node.get_leaves()
            if above:
                parent_groups.append(above)
        
        # Combine small groups if needed (for polytomies) - for display only
        child_groups = _combine_to_pairs(child_groups)
        parent_groups = _combine_to_pairs(parent_groups)
        
        if len(child_groups) < 2 or len(parent_groups) < 2:
            return None
        
        # Create quadripartition with both complete sides and sub-groups
        a, b = child_groups[0], child_groups[1]
        c, d = parent_groups[0], parent_groups[1]
        
        qp = _format_quadripartition(a, b, c, d, taxon_order, child_node, parent_node)
        
        # Add complete side information for quartet enumeration
        qp['child_side'] = child_side
        qp['parent_side'] = parent_side
        
        return qp
    
    # Traverse all parent-child relationships
    for node in root.get_all_nodes():
        if node.parent is None:
            continue
        qp = process_branch(node.parent, node)
        if qp:
            quadripartitions.append(qp)
    
    return quadripartitions


def _combine_to_pairs(groups):
    """If more than 2 groups, combine into exactly 2."""
    if len(groups) <= 2:
        return groups
    # Combine all but first into second group
    combined = groups[0]
    rest = set()
    for g in groups[1:]:
        rest.update(g)
    return [combined, rest]



# =============================================================================
# Quartet Pattern Counting
# =============================================================================

def classify_site_pattern(a, b, c, d):
    """
    Classify a site pattern for quartet (A, B, C, D).
    
    Returns:
        'AB_CD' if pattern supports (A,B)|(C,D)
        'AC_BD' if pattern supports (A,C)|(B,D)
        'AD_BC' if pattern supports (A,D)|(B,C)
        'uninformative' for constant or singleton patterns
        'missing' if any taxon has missing data
    """
    # Check for missing data
    if a in '-?' or b in '-?' or c in '-?' or d in '-?':
        return 'missing'
    
    # Check for constant patterns
    if a == b == c == d:
        return 'constant'
    
    # Check for informative 2-2 splits
    if a == b and c == d and a != c:
        return 'AB_CD'
    if a == c and b == d and a != b:
        return 'AC_BD'
    if a == d and b == c and a != b:
        return 'AD_BC'
    
    # Singleton patterns (3-1 splits) - variable but uninformative
    return 'uninformative'


def count_quartet_patterns(seq_a, seq_b, seq_c, seq_d):
    """
    Count informative site patterns supporting each quartet topology.
    
    Returns:
        (n_AB_CD, n_AC_BD, n_AD_BC)
    """
    n_AB_CD = 0
    n_AC_BD = 0
    n_AD_BC = 0
    
    nchar = len(seq_a)
    
    for i in range(nchar):
        pattern_type = classify_site_pattern(seq_a[i], seq_b[i], seq_c[i], seq_d[i])
        
        if pattern_type == 'AB_CD':
            n_AB_CD += 1
        elif pattern_type == 'AC_BD':
            n_AC_BD += 1
        elif pattern_type == 'AD_BC':
            n_AD_BC += 1
    
    return n_AB_CD, n_AC_BD, n_AD_BC



# =============================================================================
# Statistical Tests
# =============================================================================

def exact_binomial_test(successes, total, p=0.5):
    """
    Perform exact two-sided binomial test.

    Returns p-value for testing H0: probability = p.
    Uses the min-likelihood method for the two-sided test.

    All internal arithmetic is done in log space to avoid overflow when
    total is large (e.g. thousands of informative sites).
    """
    import math

    if total == 0:
        return 1.0

    log_p  = math.log(p)
    log_1p = math.log(1.0 - p)

    def log_binom_pmf(k, n):
        """log P(X = k) for Binomial(n, p), safe for large n."""
        if k < 0 or k > n:
            return -math.inf
        # log C(n,k) via lgamma to stay in floating-point
        log_comb = (math.lgamma(n + 1)
                    - math.lgamma(k + 1)
                    - math.lgamma(n - k + 1))
        return log_comb + k * log_p + (n - k) * log_1p

    log_obs = log_binom_pmf(successes, total)

    # Sum P(X = k) for all k whose log-probability <= log_obs + epsilon.
    # Accumulate in log space using log-sum-exp for numerical stability.
    log_pvalue = -math.inf
    epsilon = 1e-10  # tolerance for floating-point comparison

    for k in range(total + 1):
        lp = log_binom_pmf(k, total)
        if lp <= log_obs + epsilon:
            # log-sum-exp step
            if lp > log_pvalue:
                log_pvalue = lp + math.log1p(math.exp(log_pvalue - lp))
            else:
                log_pvalue = log_pvalue + math.log1p(math.exp(lp - log_pvalue))

    pvalue = math.exp(log_pvalue) if log_pvalue > -math.inf else 0.0
    return min(pvalue, 1.0)


def apply_multiple_test_correction(pvalues, method, alpha=0.05):
    """
    Apply multiple test correction and return adjusted p-values and significance.
    
    Methods:
    - BH: Benjamini-Hochberg (controls FDR)
    - BY: Benjamini-Yekutieli (controls FDR under arbitrary dependence)
    - bonferroni: Bonferroni correction
    - holm: Holm-Bonferroni (step-down)
    
    Returns:
        List of (adjusted_pvalue, is_significant) tuples
    """
    n = len(pvalues)
    if n == 0:
        return []
    
    # Handle NaN/None values
    valid_indices = [i for i, p in enumerate(pvalues) if p is not None and p == p]
    valid_pvalues = [pvalues[i] for i in valid_indices]
    
    if not valid_pvalues:
        return [(None, False) for _ in pvalues]
    
    m = len(valid_pvalues)
    
    # Sort by p-value
    sorted_indices = sorted(range(m), key=lambda i: valid_pvalues[i])
    sorted_pvalues = [valid_pvalues[i] for i in sorted_indices]
    
    # Calculate adjusted p-values based on method
    if method == 'bonferroni':
        # Bonferroni: multiply each p-value by total number of tests
        adjusted = [None] * m
        for i in range(m):
            adjusted[i] = min(valid_pvalues[i] * n, 1.0)
    
    elif method == 'holm':
        # Holm-Bonferroni step-down
        adjusted = [None] * m
        for rank, idx in enumerate(sorted_indices):
            adjusted[idx] = min(sorted_pvalues[rank] * (n - rank), 1.0)
        # Enforce monotonicity
        for i in range(1, m):
            adjusted[sorted_indices[i]] = max(adjusted[sorted_indices[i]], 
                                               adjusted[sorted_indices[i-1]])
    
    elif method == 'BH':
        # Benjamini-Hochberg
        adjusted = [None] * m
        prev_adj = 1.0
        for rank in range(m - 1, -1, -1):
            idx = sorted_indices[rank]
            adj = sorted_pvalues[rank] * m / (rank + 1)
            adj = min(adj, prev_adj, 1.0)
            adjusted[idx] = adj
            prev_adj = adj
    
    elif method == 'BY':
        # Benjamini-Yekutieli
        # c(m) = sum(1/i for i in 1..m)
        c_m = sum(1.0 / i for i in range(1, m + 1))
        adjusted = [None] * m
        prev_adj = 1.0
        for rank in range(m - 1, -1, -1):
            idx = sorted_indices[rank]
            adj = sorted_pvalues[rank] * m * c_m / (rank + 1)
            adj = min(adj, prev_adj, 1.0)
            adjusted[idx] = adj
            prev_adj = adj
    
    else:
        raise ValueError(f"Unknown correction method: {method}")
    
    # Map back to original indices and determine significance
    result = [(None, False) for _ in pvalues]
    for i, orig_idx in enumerate(valid_indices):
        adj_p = adjusted[i]
        result[orig_idx] = (adj_p, adj_p <= alpha)
    
    return result


# =============================================================================
# Main Analysis
# =============================================================================

def analyze_quartet_asymmetry(sequences, taxa, tree_root, taxon_order, quadripartitions):
    """
    Perform quartet asymmetry analysis.

    For each quadripartition, enumerate all VALID quartets: those with exactly
    one taxon drawn from each of the four sub-groups (a_taxa, b_taxa, c_taxa,
    d_taxa).  The tree-concordant topology is ab|cd by construction, so no
    tree-topology lookup or resolving-edge search is needed.

    Returns list of result dicts, one per valid quartet.
    """
    results = []

    for qp in quadripartitions:
        # Sort each sub-group by taxon_order for deterministic output
        a_list = sorted(qp['a_taxa'], key=lambda t: taxon_order[t])
        b_list = sorted(qp['b_taxa'], key=lambda t: taxon_order[t])
        c_list = sorted(qp['c_taxa'], key=lambda t: taxon_order[t])
        d_list = sorted(qp['d_taxa'], key=lambda t: taxon_order[t])

        node1 = qp['node1']
        node2 = qp['node2']
        qp_str = qp['formatted']

        for t_a in a_list:
            seq_a = sequences[t_a]
            for t_b in b_list:
                seq_b = sequences[t_b]
                for t_c in c_list:
                    seq_c = sequences[t_c]
                    for t_d in d_list:
                        seq_d = sequences[t_d]

                        n_ab_cd, n_ac_bd, n_ad_bc = count_quartet_patterns(
                            seq_a, seq_b, seq_c, seq_d
                        )

                        is_majority = (n_ab_cd >= n_ac_bd and n_ab_cd >= n_ad_bc)

                        minority_total = n_ac_bd + n_ad_bc
                        if minority_total > 0:
                            pvalue = exact_binomial_test(n_ac_bd, minority_total, 0.5)
                        else:
                            pvalue = 1.0

                        results.append({
                            'a': t_a,
                            'b': t_b,
                            'c': t_c,
                            'd': t_d,
                            'n_ab_cd': n_ab_cd,
                            'n_ac_bd': n_ac_bd,
                            'n_ad_bc': n_ad_bc,
                            'is_majority': is_majority,
                            'pvalue': pvalue,
                            'node1': node1,
                            'node2': node2,
                            'quadripartition': qp_str,
                        })

    return results


def summarize_by_quadripartition(results, alpha, corrections):
    """
    Summarize results by quadripartition.
    
    Returns list of summary dicts, one per unique quadripartition.
    """
    # Group results by quadripartition
    by_qp = defaultdict(list)
    for r in results:
        key = (r['node1'], r['node2'], r['quadripartition'])
        by_qp[key].append(r)
    
    summaries = []
    for (node1, node2, qp_str), qp_results in by_qp.items():
        n_quartets = len(qp_results)
        
        # Calculate means
        mean_ab_cd = sum(r['n_ab_cd'] for r in qp_results) / n_quartets
        mean_ac_bd = sum(r['n_ac_bd'] for r in qp_results) / n_quartets
        mean_ad_bc = sum(r['n_ad_bc'] for r in qp_results) / n_quartets
        
        # Count non-majority quartets
        n_non_majority = sum(1 for r in qp_results if not r['is_majority'])
        
        # Count significant at nominal alpha
        pvalues = [r['pvalue'] for r in qp_results]
        n_sig_uncorrected = sum(1 for p in pvalues if p < alpha)
        
        # Apply corrections and count significant
        sig_counts = {}
        for method in corrections:
            corrected = apply_multiple_test_correction(pvalues, method, alpha)
            n_sig = sum(1 for _, is_sig in corrected if is_sig)
            sig_counts[method] = n_sig
        
        summaries.append({
            'node1': node1,
            'node2': node2,
            'quadripartition': qp_str,
            'n_quartets': n_quartets,
            'mean_ab_cd': mean_ab_cd,
            'mean_ac_bd': mean_ac_bd,
            'mean_ad_bc': mean_ad_bc,
            'n_non_majority': n_non_majority,
            'n_sig_uncorrected': n_sig_uncorrected,
            'sig_counts': sig_counts
        })
    
    return summaries


def main():
    args = parse_args()
    
    log_progress("Quartet Asymmetry Test started")
    
    # Parse input data
    input_format = detect_format(args.input)
    log_progress("Detected input format: {}".format(input_format.upper()))
    
    if input_format == 'nexus':
        taxa, sequences = parse_nexus(args.input)
    else:
        taxa, sequences = parse_phylip(args.input)
    
    ntax = len(taxa)
    nchar = len(sequences[taxa[0]])
    log_progress("Read {} taxa, {} characters".format(ntax, nchar))
    
    # Create taxon ordering (alphabetical, case-insensitive)
    sorted_taxa = sorted(taxa, key=str.lower)
    taxon_order = {t: i for i, t in enumerate(sorted_taxa)}
    
    # Parse species tree
    log_progress("Reading species tree...")
    with open(args.tree, 'r') as f:
        tree_str = f.read().strip()
    
    tree_root = parse_newick(tree_str)
    node_map = assign_node_ids(tree_root)
    
    # Verify tree taxa match data taxa
    tree_taxa = tree_root.get_leaves()
    data_taxa = set(taxa)
    
    if tree_taxa != data_taxa:
        missing_in_tree = data_taxa - tree_taxa
        missing_in_data = tree_taxa - data_taxa
        messages = []
        if missing_in_tree:
            messages.append("  Taxa in data but not in tree: {}".format(
                ', '.join(sorted(missing_in_tree))))
        if missing_in_data:
            messages.append("  Taxa in tree but not in data: {}".format(
                ', '.join(sorted(missing_in_data))))
        print("ERROR: Taxon mismatch between data matrix and species tree:", file=sys.stderr)
        for m in messages:
            print(m, file=sys.stderr)
        print("Please correct the typo(s) and rerun.", file=sys.stderr)
        sys.exit(1)
    
    n_internal_branches = len(node_map) - 1  # Number of internal branches
    log_progress("Tree has {} internal nodes".format(len(node_map)))
    
    # Extract quadripartitions
    log_progress("Extracting quadripartitions...")
    quadripartitions = get_all_quadripartitions(tree_root, taxon_order)
    log_progress("Found {} quadripartitions".format(len(quadripartitions)))
    
    # Calculate the number of valid quartets (one taxon per sub-group per quadripartition)
    n_quartets = sum(
        len(qp['a_taxa']) * len(qp['b_taxa']) * len(qp['c_taxa']) * len(qp['d_taxa'])
        for qp in quadripartitions
    )
    log_progress("Analyzing {:,} valid quartets across {:,} quadripartitions...".format(
        n_quartets, len(quadripartitions)))
    
    # Perform analysis
    results = analyze_quartet_asymmetry(sequences, taxa, tree_root, taxon_order, quadripartitions)
    log_progress("Quartet analysis complete")
    
    # Group results by quadripartition for per-quadripartition corrections
    from collections import defaultdict
    results_by_qp = defaultdict(list)
    for i, r in enumerate(results):
        key = (r['node1'], r['node2'], r['quadripartition'])
        results_by_qp[key].append((i, r))
    
    # Apply multiple test corrections per quadripartition
    # This is more appropriate than global correction since each branch is a separate hypothesis
    corrections_applied = {method: [None] * len(results) for method in args.correction}
    
    for qp_key, qp_results in results_by_qp.items():
        # Extract indices and p-values for this quadripartition
        indices = [i for i, r in qp_results]
        pvalues = [r['pvalue'] for i, r in qp_results]
        
        # Apply corrections within this quadripartition
        for method in args.correction:
            corrected = apply_multiple_test_correction(pvalues, method, args.alpha)
            for j, orig_idx in enumerate(indices):
                corrections_applied[method][orig_idx] = corrected[j]
    
    # Write detailed results
    detail_file = "{}.quartet_asym.tsv".format(args.output)
    log_progress("Writing detailed results to {}".format(detail_file))
    
    with open(detail_file, 'w') as f:
        # Header
        header = ['taxon_a', 'taxon_b', 'taxon_c', 'taxon_d', 
                  'n_ab_cd', 'n_ac_bd', 'n_ad_bc', 'is_majority', 'pvalue']
        for method in args.correction:
            header.extend(['adj_p_{}'.format(method), 'sig_{}'.format(method)])
        header.extend(['node1', 'node2', 'quadripartition'])
        f.write('\t'.join(header) + '\n')
        
        for i, r in enumerate(results):
            row = [
                r['a'], r['b'], r['c'], r['d'],
                str(r['n_ab_cd']), str(r['n_ac_bd']), str(r['n_ad_bc']),
                'yes' if r['is_majority'] else 'no',
                '{:.6g}'.format(r['pvalue'])
            ]
            for method in args.correction:
                adj_p, is_sig = corrections_applied[method][i]
                row.append('{:.6g}'.format(adj_p) if adj_p is not None else 'NA')
                row.append('yes' if is_sig else 'no')
            row.extend([str(r['node1']), str(r['node2']), r['quadripartition']])
            f.write('\t'.join(row) + '\n')
    
    # Write summary by quadripartition
    summaries = summarize_by_quadripartition(results, args.alpha, args.correction)
    summary_file = "{}.quadripartition_summary.tsv".format(args.output)
    log_progress("Writing summary to {}".format(summary_file))
    
    with open(summary_file, 'w') as f:
        header = ['node1', 'node2', 'n_quartets', 
                  'mean_ab_cd', 'mean_ac_bd', 'mean_ad_bc',
                  'n_non_majority', 'n_sig_uncorrected']
        for method in args.correction:
            header.append('n_sig_{}'.format(method))
        header.append('quadripartition')
        f.write('\t'.join(header) + '\n')
        
        for s in summaries:
            row = [
                str(s['node1']), str(s['node2']), str(s['n_quartets']),
                '{:.2f}'.format(s['mean_ab_cd']),
                '{:.2f}'.format(s['mean_ac_bd']),
                '{:.2f}'.format(s['mean_ad_bc']),
                str(s['n_non_majority']),
                str(s['n_sig_uncorrected'])
            ]
            for method in args.correction:
                row.append(str(s['sig_counts'].get(method, 0)))
            row.append(s['quadripartition'])
            f.write('\t'.join(row) + '\n')
    
    # Write labeled tree
    tree_file = "{}.labeled_tree.tre".format(args.output)
    log_progress("Writing labeled tree to {}".format(tree_file))
    
    with open(tree_file, 'w') as f:
        newick_str = tree_to_newick(tree_root, include_node_ids=True, 
                                     include_branch_lengths=True)
        f.write(newick_str + ';\n')
    
    # Print summary statistics
    print("\n" + "="*70)
    print("Quartet Asymmetry Test Summary")
    print("="*70)
    print("Input file:        {}".format(args.input))
    print("Species tree:      {}".format(args.tree))
    print("Taxa:              {}".format(ntax))
    print("Characters:        {}".format(nchar))
    print("Quadripartitions:  {}".format(len(quadripartitions)))
    print("Total quartets:    {}".format(len(results)))
    print("Alpha level:       {}".format(args.alpha))
    print("-"*70)
    
    # Count overall statistics
    n_majority = sum(1 for r in results if r['is_majority'])
    n_non_majority = len(results) - n_majority
    n_sig_uncorrected = sum(1 for r in results if r['pvalue'] < args.alpha)
    
    print("Majority quartets (ab|cd highest): {} ({:.1f}%)".format(
        n_majority, 100 * n_majority / len(results)))
    print("Non-majority quartets:             {} ({:.1f}%)".format(
        n_non_majority, 100 * n_non_majority / len(results)))
    print("-"*70)
    print("Asymmetry tests (minority quartets):")
    print("  Significant uncorrected (P<{}): {} ({:.1f}%)".format(
        args.alpha, n_sig_uncorrected, 100 * n_sig_uncorrected / len(results)))

    for method in args.correction:
        n_sig = sum(1 for _, is_sig in corrections_applied[method] if is_sig)
        print("  Significant {} corrected:    {} ({:.1f}%)".format(
            method.ljust(8), n_sig, 100 * n_sig / len(results)))

    print("-"*70)
    print("Branches with at least one significant minority quartet asymmetry:")

    n_branches = len(quadripartitions)

    # Uncorrected: branches with at least one p < alpha
    branches_uncorrected = set(
        (r['node1'], r['node2']) for r in results if r['pvalue'] < args.alpha
    )
    print("  Significant uncorrected (P<{}): {} of {} branches".format(
        args.alpha, len(branches_uncorrected), n_branches))

    # Corrected: branches with at least one significant corrected result
    for method in args.correction:
        branches_sig = set(
            (results[i]['node1'], results[i]['node2'])
            for i, (_, is_sig) in enumerate(corrections_applied[method])
            if is_sig
        )
        print("  Significant {} corrected:    {} of {} branches".format(
            method.ljust(8), len(branches_sig), n_branches))

    print("-"*70)
    print("Output files:")
    print("  Detailed results: {}".format(detail_file))
    print("  Summary:          {}".format(summary_file))
    print("  Labeled tree:     {}".format(tree_file))
    print("="*70)
    
    log_progress("Analysis complete")


if __name__ == '__main__':
    main()
