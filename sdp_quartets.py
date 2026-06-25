#!/usr/bin/env python3
"""
sdp_quartets.py - Split Decomposition with Parsimony for Quartets

A quartet-based phylogenetic analysis tool that:
1. Enumerates all C(n,4) quartets from a binary character matrix
2. Determines the most parsimonious quartet tree(s) by pattern counting
3. Generates an MRP (Matrix Representation using Parsimony) matrix
4. Outputs a detailed quartet log with scores and resolution status
5. Optionally performs bootstrap analysis with weighted consensus

Based on the method described in Springer et al. 2020 (doi:10.1093/jhered/esz076)
"""

import argparse
import sys
import re
import os
import shutil
import subprocess
import random
import time
from datetime import datetime
from itertools import combinations


def timestamp():
    """Return current timestamp string."""
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def log_progress(message):
    """Print timestamped progress message to stderr."""
    print("[{}] {}".format(timestamp(), message), file=sys.stderr)


def parse_args():
    parser = argparse.ArgumentParser(
        description='SDPquartets: Split Decomposition with Parsimony for Quartets',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python sdp_quartets.py --input data.phy --output results
  python sdp_quartets.py --input data.nex --output results --format phylip
  python sdp_quartets.py --input data.phy --output results --paup
  python sdp_quartets.py --input data.phy --output results --paup --paup-path /usr/local/bin/paup
  python sdp_quartets.py --input data.phy --output results --springer-compat
  python sdp_quartets.py --input data.phy --output results --bootstrap --reps 100 --paup
  python sdp_quartets.py --input data.phy --output results -b --reps 500 --seed 12345 --paup
  python sdp_quartets.py --input data.phy --output results -b --boot-no-opt --lazy --reps 100 --paup
  python sdp_quartets.py --input data.phy --output results -b --reps 100 --clean --paup
        """
    )
    parser.add_argument('--input', '-i', required=True,
                        help='Input character matrix (PHYLIP or NEXUS format)')
    parser.add_argument('--output', '-o', required=True,
                        help='Base name for output files')
    parser.add_argument('--format', '-f', choices=['phylip', 'nexus'], default='nexus',
                        help='Output MRP matrix format (default: nexus)')
    parser.add_argument('--springer-compat', action='store_true',
                        help='Include 3-way ties with 6x/3x/2x weighting for compatibility '
                             'with Springer et al. 2020')
    parser.add_argument('--paup', action='store_true',
                        help='Run PAUP* tree search on MRP matrix')
    parser.add_argument('--paup-path', default='paup',
                        help='Path to PAUP* executable (default: paup)')
    parser.add_argument('--bootstrap', '-b', action='store_true',
                        help='Perform bootstrap analysis')
    parser.add_argument('--reps', type=int, default=100,
                        help='Number of bootstrap replicates (default: 100)')
    parser.add_argument('--seed', type=int, default=None,
                        help='Random number seed for bootstrap (default: from system clock)')
    parser.add_argument('--boot-no-opt', action='store_true',
                        help='Skip optimal tree search when bootstrapping (only do bootstrap)')
    parser.add_argument('--lazy', action='store_true',
                        help='Use faster but less thorough PAUP* search settings')
    parser.add_argument('--clean', action='store_true',
                        help='Remove bootstrap directory and intermediate files after analysis')
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


def bootstrap_resample(sequences, taxa):
    """
    Create a bootstrap-resampled version of the sequence data.
    Sample characters with replacement.
    """
    nchar = len(sequences[taxa[0]])
    # Sample indices with replacement
    indices = [random.randint(0, nchar - 1) for _ in range(nchar)]
    
    # Create new sequences
    resampled = {}
    for taxon in taxa:
        orig_seq = sequences[taxon]
        resampled[taxon] = ''.join(orig_seq[i] for i in indices)
    
    return resampled


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
    
    pattern = (a, b, c, d)
    
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


def score_quartet(seq_a, seq_b, seq_c, seq_d):
    """
    Calculate parsimony scores for all three quartet topologies.
    
    Returns:
        (score_AB_CD, score_AC_BD, score_AD_BC, n_informative, n_uninformative)
    """
    n_AB_CD = 0  # patterns supporting (A,B)|(C,D)
    n_AC_BD = 0  # patterns supporting (A,C)|(B,D)
    n_AD_BC = 0  # patterns supporting (A,D)|(B,C)
    n_uninformative = 0  # variable but uninformative (singletons, missing)
    
    nchar = len(seq_a)
    
    for i in range(nchar):
        pattern_type = classify_site_pattern(seq_a[i], seq_b[i], seq_c[i], seq_d[i])
        
        if pattern_type == 'AB_CD':
            n_AB_CD += 1
        elif pattern_type == 'AC_BD':
            n_AC_BD += 1
        elif pattern_type == 'AD_BC':
            n_AD_BC += 1
        elif pattern_type == 'uninformative' or pattern_type == 'missing':
            n_uninformative += 1
        # 'constant' patterns contribute 0 steps, ignore
    
    # Calculate scores
    # Each informative pattern: 1 step for concordant topology, 2 for discordant
    # Uninformative variable patterns: 1 step for all topologies
    total_informative = n_AB_CD + n_AC_BD + n_AD_BC
    
    score_AB_CD = n_AB_CD + 2 * (n_AC_BD + n_AD_BC) + n_uninformative
    score_AC_BD = n_AC_BD + 2 * (n_AB_CD + n_AD_BC) + n_uninformative
    score_AD_BC = n_AD_BC + 2 * (n_AB_CD + n_AC_BD) + n_uninformative
    
    return score_AB_CD, score_AC_BD, score_AD_BC, total_informative, n_uninformative


def determine_resolution(score_AB_CD, score_AC_BD, score_AD_BC):
    """
    Determine the resolution status and winning topology/topologies.
    
    Returns:
        (resolution, winning_topologies)
        resolution: 1 (resolved), 2 (two-way tie), 3 (three-way tie)
        winning_topologies: list of topology strings
    """
    min_score = min(score_AB_CD, score_AC_BD, score_AD_BC)
    
    winners = []
    if score_AB_CD == min_score:
        winners.append('AB_CD')
    if score_AC_BD == min_score:
        winners.append('AC_BD')
    if score_AD_BC == min_score:
        winners.append('AD_BC')
    
    return len(winners), winners


def topology_to_newick(topology, taxa):
    """Convert topology code to Newick string."""
    a, b, c, d = taxa
    if topology == 'AB_CD':
        return "(({},{}),({},{}))".format(a, b, c, d)
    elif topology == 'AC_BD':
        return "(({},{}),({},{}))".format(a, c, b, d)
    elif topology == 'AD_BC':
        return "(({},{}),({},{}))".format(a, d, b, c)


def generate_mrp_character(topology, quartet_taxa, all_taxa):
    """
    Generate MRP character for a quartet topology.
    
    Returns a dictionary mapping taxon -> character state ('0', '1', or '?')
    """
    a, b, c, d = quartet_taxa
    char = {taxon: '?' for taxon in all_taxa}
    
    if topology == 'AB_CD':
        char[a] = '1'
        char[b] = '1'
        char[c] = '0'
        char[d] = '0'
    elif topology == 'AC_BD':
        char[a] = '1'
        char[c] = '1'
        char[b] = '0'
        char[d] = '0'
    elif topology == 'AD_BC':
        char[a] = '1'
        char[d] = '1'
        char[b] = '0'
        char[c] = '0'
    
    return char


def write_mrp_phylip(filename, taxa, mrp_chars):
    """Write MRP matrix in relaxed PHYLIP format."""
    if not mrp_chars:
        log_progress("Warning: No MRP characters to write")
        return
    
    nchar = len(mrp_chars)
    ntax = len(taxa)
    
    with open(filename, 'w') as f:
        f.write("{} {}\n".format(ntax, nchar))
        for taxon in taxa:
            seq = ''.join(char[taxon] for char in mrp_chars)
            f.write("{}\t{}\n".format(taxon, seq))


def write_mrp_nexus(filename, taxa, mrp_chars):
    """Write MRP matrix in NEXUS format."""
    if not mrp_chars:
        log_progress("Warning: No MRP characters to write")
        return
    
    nchar = len(mrp_chars)
    ntax = len(taxa)
    
    with open(filename, 'w') as f:
        f.write("#NEXUS\n")
        f.write("Begin data;\n")
        f.write("Dimensions ntax={} nchar={};\n".format(ntax, nchar))
        f.write("Format datatype=standard symbols=\"01\" missing=?;\n")
        f.write("Matrix\n")
        for taxon in taxa:
            seq = ''.join(char[taxon] for char in mrp_chars)
            f.write("{}\t{}\n".format(taxon, seq))
        f.write(";\nEnd;\n")


def check_paup_available(paup_path):
    """Check if PAUP* is available at the specified path."""
    try:
        result = subprocess.run(
            [paup_path, '-h'],
            capture_output=True,
            text=True,
            timeout=10
        )
        # PAUP -h returns help text; we just need it to run without error
        return True
    except FileNotFoundError:
        return False
    except subprocess.TimeoutExpired:
        return False
    except Exception as e:
        log_progress("Warning: Error checking PAUP* availability: {}".format(e))
        return False


def append_paup_block(filename, output_base, include_contree=True, lazy=False, rseed=None):
    """
    Append PAUP block to NEXUS file for tree search.
    
    Args:
        filename: NEXUS file to append to
        output_base: Base name for output tree files
        include_contree: Whether to include consensus tree command
        lazy: Use faster but less thorough search settings
        rseed: Random seed for PAUP (0-999999999)
    """
    if rseed is None:
        rseed = random.randint(0, 999999999)
    
    if lazy:
        # Fast search: limited maxtrees, chuck suboptimal trees
        if include_contree:
            paup_block = """
Begin PAUP;
    set maxtrees = 1000 increase=no;
    hsearch addseq=random nreps=100 nchuck=10 chuckscore=1 rseed={};
    savetrees file={}.SDPtrees.tre format=newick trees=all replace;
    contree all / treefile={}.SDPconstree.tre replace;
    quit;
End;
""".format(rseed, output_base, output_base)
        else:
            paup_block = """
Begin PAUP;
    set maxtrees = 1000 increase=no;
    hsearch addseq=random nreps=100 nchuck=10 chuckscore=1 rseed={};
    savetrees file={}.SDPtrees.tre format=newick trees=all replace;
    quit;
End;
""".format(rseed, output_base)
    else:
        # Thorough search: autoincrement maxtrees
        if include_contree:
            paup_block = """
Begin PAUP;
    set increase=auto;
    hsearch addseq=random nreps=100 rseed={};
    savetrees file={}.SDPtrees.tre format=newick trees=all replace;
    contree all / treefile={}.SDPconstree.tre replace;
    quit;
End;
""".format(rseed, output_base, output_base)
        else:
            paup_block = """
Begin PAUP;
    set increase=auto;
    hsearch addseq=random nreps=100 rseed={};
    savetrees file={}.SDPtrees.tre format=newick trees=all replace;
    quit;
End;
""".format(rseed, output_base)
    
    with open(filename, 'a') as f:
        f.write(paup_block)


def run_paup(paup_path, nexus_file, quiet=False, cwd=None):
    """Execute PAUP* on the NEXUS file."""
    try:
        result = subprocess.run(
            [paup_path, '-n', nexus_file],
            capture_output=True,
            text=True,
            cwd=cwd
        )
        if result.returncode != 0 and not quiet:
            log_progress("Warning: PAUP* returned non-zero exit code: {}".format(result.returncode))
            if result.stderr:
                log_progress("PAUP* stderr: {}".format(result.stderr))
        return result.returncode == 0
    except Exception as e:
        if not quiet:
            log_progress("Error running PAUP*: {}".format(e))
        return False


def read_newick_trees(filename):
    """Read Newick trees from a file, one per line."""
    trees = []
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if line and line.startswith('('):
                    # Remove trailing semicolon if present
                    if line.endswith(';'):
                        line = line[:-1]
                    trees.append(line)
    except FileNotFoundError:
        pass
    return trees


def analyze_quartets(sequences, taxa, springer_compat, write_log=True, log_file=None):
    """
    Analyze all quartets and generate MRP characters.
    
    Returns:
        (mrp_chars, quartet_trees, count_resolved, count_2way, count_3way)
    """
    mrp_chars = []
    quartet_trees = []
    count_resolved = 0
    count_2way = 0
    count_3way = 0
    
    sorted_taxa = sorted(taxa, key=str.lower)
    
    for quartet in combinations(sorted_taxa, 4):
        a, b, c, d = quartet
        
        # Get scores
        score_AB_CD, score_AC_BD, score_AD_BC, n_inf, n_uninf = score_quartet(
            sequences[a], sequences[b], sequences[c], sequences[d]
        )
        
        # Determine resolution
        resolution, winners = determine_resolution(score_AB_CD, score_AC_BD, score_AD_BC)
        
        # Write to log if requested
        if write_log and log_file:
            log_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                a, b, c, d, resolution, score_AB_CD, score_AC_BD, score_AD_BC))
        
        # Update statistics
        if resolution == 1:
            count_resolved += 1
        elif resolution == 2:
            count_2way += 1
        else:
            count_3way += 1
        
        # Generate trees and MRP characters based on mode
        if springer_compat:
            # Springer et al. weighting: 6x/3x/2x per topology
            # Resolved: 6 chars, 2-way: 3+3 chars, 3-way: 2+2+2 chars
            if resolution == 1:
                weight = 6
            elif resolution == 2:
                weight = 3
            else:
                weight = 2
            
            for topo in winners:
                newick = topology_to_newick(topo, quartet)
                for _ in range(weight):
                    quartet_trees.append(newick)
                    mrp_chars.append(generate_mrp_character(topo, quartet, taxa))
        else:
            # Default mode: compact weighting that maintains 2:1 ratio
            # Resolved: 2 chars, 2-way: 1+1 chars, 3-way: excluded
            if resolution == 1:
                # Single winner - weight 2
                topo = winners[0]
                newick = topology_to_newick(topo, quartet)
                for _ in range(2):
                    quartet_trees.append(newick)
                    mrp_chars.append(generate_mrp_character(topo, quartet, taxa))
            elif resolution == 2:
                # Two-way tie - weight 1 each (informative by exclusion)
                for topo in winners:
                    newick = topology_to_newick(topo, quartet)
                    quartet_trees.append(newick)
                    mrp_chars.append(generate_mrp_character(topo, quartet, taxa))
            # 3-way ties: skip (uninformative)
    
    return mrp_chars, quartet_trees, count_resolved, count_2way, count_3way


def run_bootstrap_replicate(rep_num, sequences, taxa, springer_compat, output_dir, 
                            output_base, paup_path, lazy=False):
    """
    Run a single bootstrap replicate.
    
    Returns:
        (trees_with_weights, n_trees) where trees_with_weights is a list of 
        (tree, weight) tuples and n_trees is the count of MP trees found
    """
    # Bootstrap resample
    resampled = bootstrap_resample(sequences, taxa)
    
    # Analyze quartets (no log file)
    mrp_chars, quartet_trees, _, _, _ = analyze_quartets(
        resampled, taxa, springer_compat, write_log=False
    )
    
    # Write MRP matrix - use full path
    rep_str = "rep{:03d}".format(rep_num)
    mrp_basename = "{}.{}.mrp.nex".format(output_base, rep_str)
    mrp_filename = os.path.join(output_dir, mrp_basename)
    write_mrp_nexus(mrp_filename, taxa, mrp_chars)
    
    # Generate random seed for this replicate's PAUP search
    rseed = random.randint(0, 999999999)
    
    # Append PAUP block - use basename only (PAUP will run from output_dir)
    rep_basename_noext = "{}.{}".format(output_base, rep_str)
    append_paup_block(mrp_filename, rep_basename_noext, include_contree=False, 
                      lazy=lazy, rseed=rseed)
    
    # Run PAUP from the bootstrap directory so output files go there
    success = run_paup(paup_path, mrp_basename, quiet=True, cwd=output_dir)
    
    if not success:
        log_progress("Warning: PAUP* failed for replicate {}".format(rep_num))
        return [], 0
    
    # Read resulting trees from the bootstrap directory
    tree_filename = os.path.join(output_dir, "{}.SDPtrees.tre".format(rep_basename_noext))
    trees = read_newick_trees(tree_filename)
    
    if not trees:
        log_progress("Warning: No trees found for replicate {}".format(rep_num))
        return [], 0
    
    # Calculate weight (1/n for n trees)
    n_trees = len(trees)
    if n_trees > 1:
        weight = "1/{}".format(n_trees)
    else:
        weight = "1"
    
    return [(tree, weight) for tree in trees], n_trees


def write_weighted_treefile(filename, all_trees_by_rep):
    """
    Write NEXUS treefile with weighted trees for bootstrap consensus.
    
    all_trees_by_rep: list of lists, where each inner list contains (tree, weight) tuples
    """
    with open(filename, 'w') as f:
        f.write("#NEXUS\n")
        f.write("Begin PAUP; set increase=auto; End;\n")
        f.write("Begin trees;\n")
        
        for rep_num, trees in enumerate(all_trees_by_rep, 1):
            if not trees:
                continue
            f.write("[Trees found in bootstrap replicate #{}]\n".format(rep_num))
            for tree_idx, (tree, weight) in enumerate(trees, 1):
                f.write("tree B_{}.{} = [&W {}] {};\n".format(rep_num, tree_idx, weight, tree))
        
        f.write("End;\n")


def compute_bootstrap_stats(trees_per_rep):
    """
    Compute statistics on the number of MP trees per bootstrap replicate.
    
    Args:
        trees_per_rep: list of tree counts for each replicate
    
    Returns:
        dict with keys: single_tree_reps, multi_tree_reps, min, q1, median, q3, max
    """
    # Filter out failed replicates (0 trees)
    valid_counts = [n for n in trees_per_rep if n > 0]
    
    if not valid_counts:
        return {
            'single_tree_reps': 0,
            'multi_tree_reps': 0,
            'min': 0,
            'q1': 0,
            'median': 0,
            'q3': 0,
            'max': 0
        }
    
    # Count single vs multiple tree replicates
    single_tree_reps = sum(1 for n in valid_counts if n == 1)
    multi_tree_reps = sum(1 for n in valid_counts if n > 1)
    
    # Sort for percentile calculations
    sorted_counts = sorted(valid_counts)
    n = len(sorted_counts)
    
    # Min and max
    min_val = sorted_counts[0]
    max_val = sorted_counts[-1]
    
    # Median (Q2)
    if n % 2 == 0:
        median = (sorted_counts[n // 2 - 1] + sorted_counts[n // 2]) / 2.0
    else:
        median = sorted_counts[n // 2]
    
    # Q1 (25th percentile) - use linear interpolation method
    q1_pos = (n - 1) * 0.25
    q1_lower = int(q1_pos)
    q1_frac = q1_pos - q1_lower
    if q1_lower + 1 < n:
        q1 = sorted_counts[q1_lower] * (1 - q1_frac) + sorted_counts[q1_lower + 1] * q1_frac
    else:
        q1 = sorted_counts[q1_lower]
    
    # Q3 (75th percentile)
    q3_pos = (n - 1) * 0.75
    q3_lower = int(q3_pos)
    q3_frac = q3_pos - q3_lower
    if q3_lower + 1 < n:
        q3 = sorted_counts[q3_lower] * (1 - q3_frac) + sorted_counts[q3_lower + 1] * q3_frac
    else:
        q3 = sorted_counts[q3_lower]
    
    return {
        'single_tree_reps': single_tree_reps,
        'multi_tree_reps': multi_tree_reps,
        'min': min_val,
        'q1': q1,
        'median': median,
        'q3': q3,
        'max': max_val
    }


def main():
    args = parse_args()
    
    log_progress("SDPquartets analysis started")
    
    # Handle random seed
    seed = None
    if args.seed is not None:
        seed = args.seed
        log_progress("Using provided random seed: {}".format(seed))
        random.seed(seed)
    elif args.bootstrap or args.paup:
        # Generate seed for reproducibility
        seed = int(time.time() * 1000) % (2**31)
        log_progress("Generated random seed from system clock: {}".format(seed))
        random.seed(seed)
    
    # Check PAUP availability if --paup or --bootstrap specified
    if args.paup or args.bootstrap:
        if not check_paup_available(args.paup_path):
            log_progress("Error: PAUP* not found at '{}'".format(args.paup_path))
            log_progress("Please install PAUP* or specify the correct path with --paup-path")
            sys.exit(1)
        log_progress("PAUP* found at: {}".format(args.paup_path))
        
        # Force NEXUS format if running PAUP
        if args.format != 'nexus':
            log_progress("Note: Switching to NEXUS format for PAUP* analysis")
            args.format = 'nexus'
    
    # Bootstrap requires PAUP
    if args.bootstrap and not args.paup:
        log_progress("Note: Bootstrap analysis requires PAUP*, enabling --paup")
        args.paup = True
    
    # Detect and parse input format
    input_format = detect_format(args.input)
    log_progress("Detected input format: {}".format(input_format.upper()))
    
    if input_format == 'nexus':
        taxa, sequences = parse_nexus(args.input)
    else:
        taxa, sequences = parse_phylip(args.input)
    
    ntax = len(taxa)
    nchar = len(sequences[taxa[0]])
    log_progress("Read {} taxa, {} characters".format(ntax, nchar))
    
    # Calculate number of quartets
    n_quartets = ntax * (ntax - 1) * (ntax - 2) * (ntax - 3) // 24
    log_progress("Analyzing {} quartets...".format(n_quartets))
    
    # Open quartet log file (only for main analysis, not bootstrap)
    log_filename = "{}.quartet_log.tsv".format(args.output)
    log_file = open(log_filename, 'w')
    log_file.write("TaxonA\tTaxonB\tTaxonC\tTaxonD\tResolution\tScore_AB_CD\tScore_AC_BD\tScore_AD_BC\n")
    
    # Analyze quartets
    mrp_chars, quartet_trees, count_resolved, count_2way, count_3way = analyze_quartets(
        sequences, taxa, args.springer_compat, write_log=True, log_file=log_file
    )
    
    log_file.close()
    log_progress("Quartet analysis complete: {} resolved, {} 2-way ties, {} 3-way ties".format(
        count_resolved, count_2way, count_3way))
    
    # Write quartet trees
    tree_filename = "{}.quartets.tre".format(args.output)
    with open(tree_filename, 'w') as f:
        for tree in quartet_trees:
            f.write(tree + ";\n")
    
    # Write MRP matrix
    if args.format == 'nexus':
        mrp_filename = "{}.mrp.nex".format(args.output)
        write_mrp_nexus(mrp_filename, taxa, mrp_chars)
    else:
        mrp_filename = "{}.mrp.phy".format(args.output)
        write_mrp_phylip(mrp_filename, taxa, mrp_chars)
    
    log_progress("MRP matrix written: {} characters".format(len(mrp_chars)))
    
    # Run PAUP if requested (unless --boot-no-opt is set with bootstrap)
    paup_success = False
    skip_optimal = args.bootstrap and args.boot_no_opt
    
    if args.paup and not args.bootstrap:
        # Standard PAUP search (no bootstrap)
        log_progress("Appending PAUP* block and running tree search...")
        rseed = random.randint(0, 999999999)
        append_paup_block(mrp_filename, args.output, include_contree=True, 
                          lazy=args.lazy, rseed=rseed)
        paup_success = run_paup(args.paup_path, mrp_filename)
        if paup_success:
            log_progress("PAUP* analysis completed successfully")
        else:
            log_progress("Warning: PAUP* analysis may have encountered issues")
    elif args.paup and args.bootstrap and not skip_optimal:
        # Bootstrap mode but still want optimal tree search first
        log_progress("Appending PAUP* block and running optimal tree search...")
        rseed = random.randint(0, 999999999)
        append_paup_block(mrp_filename, args.output, include_contree=True,
                          lazy=args.lazy, rseed=rseed)
        paup_success = run_paup(args.paup_path, mrp_filename)
        if paup_success:
            log_progress("PAUP* optimal tree analysis completed successfully")
        else:
            log_progress("Warning: PAUP* optimal tree analysis may have encountered issues")
    
    # Bootstrap analysis
    bootstrap_dir = None
    bootstrap_stats = None
    if args.bootstrap:
        if args.lazy:
            log_progress("Starting bootstrap analysis with {} replicates (lazy search mode)...".format(args.reps))
        else:
            log_progress("Starting bootstrap analysis with {} replicates...".format(args.reps))
        
        # Create bootstrap directory
        bootstrap_dir = args.output
        os.makedirs(bootstrap_dir, exist_ok=True)
        log_progress("Bootstrap files will be written to: {}/".format(bootstrap_dir))
        
        all_trees_by_rep = []
        trees_per_rep = []
        
        for rep in range(1, args.reps + 1):
            if rep % 10 == 0 or rep == 1 or rep == args.reps:
                log_progress("Processing bootstrap replicate {}/{}".format(rep, args.reps))
            
            trees, n_trees = run_bootstrap_replicate(
                rep, sequences, taxa, args.springer_compat,
                bootstrap_dir, args.output, args.paup_path, lazy=args.lazy
            )
            all_trees_by_rep.append(trees)
            trees_per_rep.append(n_trees)
        
        # Write combined weighted treefile
        weighted_treefile = "{}.bootstrap_trees.tre".format(args.output)
        write_weighted_treefile(weighted_treefile, all_trees_by_rep)
        log_progress("Bootstrap weighted treefile written: {}".format(weighted_treefile))
        
        # Compute bootstrap statistics
        bootstrap_stats = compute_bootstrap_stats(trees_per_rep)
        
        # Count total trees
        total_trees = sum(len(trees) for trees in all_trees_by_rep)
        successful_reps = sum(1 for trees in all_trees_by_rep if trees)
        log_progress("Bootstrap complete: {}/{} successful replicates, {} total trees".format(
            successful_reps, args.reps, total_trees))
        
        # Clean up bootstrap directory if requested
        if args.clean:
            try:
                shutil.rmtree(bootstrap_dir)
                log_progress("Cleaned up bootstrap directory: {}/".format(bootstrap_dir))
            except Exception as e:
                log_progress("Warning: Could not remove bootstrap directory: {}".format(e))
    
    # Summary statistics
    print("\n" + "="*60, file=sys.stderr)
    print("SDPquartets Analysis Summary", file=sys.stderr)
    print("="*60, file=sys.stderr)
    print("Input file:        {}".format(args.input), file=sys.stderr)
    print("Taxa:              {}".format(ntax), file=sys.stderr)
    print("Characters:        {}".format(nchar), file=sys.stderr)
    print("Total quartets:    {}".format(n_quartets), file=sys.stderr)
    print("-"*60, file=sys.stderr)
    print("Resolved (1 MP):   {} ({:.2f}%)".format(count_resolved, 100*count_resolved/n_quartets), 
          file=sys.stderr)
    print("Two-way ties:      {} ({:.2f}%)".format(count_2way, 100*count_2way/n_quartets), 
          file=sys.stderr)
    print("Three-way ties:    {} ({:.2f}%)".format(count_3way, 100*count_3way/n_quartets), 
          file=sys.stderr)
    print("-"*60, file=sys.stderr)
    print("MRP characters:    {}".format(len(mrp_chars)), file=sys.stderr)
    if args.springer_compat:
        print("Mode:              Springer et al. compatible (6x/3x/2x weighting)", 
              file=sys.stderr)
    else:
        print("Mode:              Compact (2x/1x, 3-way ties excluded)", file=sys.stderr)
    if args.lazy:
        print("Search mode:       Lazy (fast, may miss some optimal trees)", file=sys.stderr)
    else:
        print("Search mode:       Thorough (increase=yes)", file=sys.stderr)
    print("-"*60, file=sys.stderr)
    print("Output files:", file=sys.stderr)
    print("  Quartet log:     {}".format(log_filename), file=sys.stderr)
    print("  Quartet trees:   {}".format(tree_filename), file=sys.stderr)
    print("  MRP matrix:      {}".format(mrp_filename), file=sys.stderr)
    if args.paup and paup_success and not args.bootstrap:
        print("  PAUP* trees:     {}.SDPtrees.tre".format(args.output), file=sys.stderr)
        print("  Consensus tree:  {}.SDPconstree.tre".format(args.output), file=sys.stderr)
        if seed is not None:
            print("  Random seed:     {}".format(seed), file=sys.stderr)
    if args.bootstrap:
        if paup_success and not skip_optimal:
            print("  PAUP* trees:     {}.SDPtrees.tre".format(args.output), file=sys.stderr)
            print("  Consensus tree:  {}.SDPconstree.tre".format(args.output), file=sys.stderr)
        if args.clean:
            print("  Bootstrap dir:   {}/ (cleaned)".format(bootstrap_dir), file=sys.stderr)
        else:
            print("  Bootstrap dir:   {}/".format(bootstrap_dir), file=sys.stderr)
        print("  Bootstrap trees: {}.bootstrap_trees.tre".format(args.output), file=sys.stderr)
        print("  Random seed:     {}".format(seed), file=sys.stderr)
        print("-"*60, file=sys.stderr)
        print("Bootstrap statistics (MP trees per replicate):", file=sys.stderr)
        print("  Single-tree replicates:   {}".format(bootstrap_stats['single_tree_reps']), 
              file=sys.stderr)
        print("  Multi-tree replicates:    {}".format(bootstrap_stats['multi_tree_reps']), 
              file=sys.stderr)
        # Format quartile values - show as integers if whole numbers, else 1 decimal
        def fmt_stat(val):
            if val == int(val):
                return str(int(val))
            else:
                return "{:.1f}".format(val)
        print("  Min / Q1 / Median / Q3 / Max: {} / {} / {} / {} / {}".format(
            fmt_stat(bootstrap_stats['min']),
            fmt_stat(bootstrap_stats['q1']),
            fmt_stat(bootstrap_stats['median']),
            fmt_stat(bootstrap_stats['q3']),
            fmt_stat(bootstrap_stats['max'])), file=sys.stderr)
    print("="*60, file=sys.stderr)
    log_progress("Analysis complete")


if __name__ == '__main__':
    main()
