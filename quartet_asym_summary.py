#!/usr/bin/env python3
"""
quartet_asym_summary.py - Taxon involvement summary for quartet asymmetry results

Reads the per-quartet output from rgc_quartet_asym.py and produces a text
summary of taxon involvement in significant asymmetric quartets.

For each branch with at least one significant quartet, reports:
  - The quadripartition
  - Direction of the asymmetry (ac|bd vs ad|bc favored)
  - Binomial test for directional bias
  - Per-taxon counts and bias scores

Output:
  OUTPREFIX.summary.txt

No external dependencies beyond the Python standard library.

Usage:
  python quartet_asym_summary.py -i results.quartet_asym.tsv -o OUTPREFIX
  python quartet_asym_summary.py -i results.quartet_asym.tsv -o OUTPREFIX \\
      --correction BH --alpha 0.05
"""

import argparse
import csv
import math
import re
import sys
from collections import Counter, defaultdict


# =============================================================================
# Quadripartition parsing
# =============================================================================

def parse_quadripartition(qp_str):
    """
    Parse a quadripartition string into its four taxon groups.

    Format: (group_a),(group_b)|(group_c),(group_d)

    Returns (taxon_to_group, a_taxa, b_taxa, c_taxa, d_taxa).
    """
    if not qp_str or qp_str == 'unresolved':
        return {}, set(), set(), set(), set()

    sides = qp_str.split('|')
    if len(sides) != 2:
        return {}, set(), set(), set(), set()

    def parse_side(s):
        groups = re.findall(r'\(([^)]+)\)', s)
        return groups if len(groups) == 2 else ['', '']

    left_groups  = parse_side(sides[0])
    right_groups = parse_side(sides[1])

    a_taxa = set(t.strip() for t in left_groups[0].split(',')  if t.strip())
    b_taxa = set(t.strip() for t in left_groups[1].split(',')  if t.strip())
    c_taxa = set(t.strip() for t in right_groups[0].split(',') if t.strip())
    d_taxa = set(t.strip() for t in right_groups[1].split(',') if t.strip())

    taxon_to_group = {}
    for t in a_taxa: taxon_to_group[t] = 'a'
    for t in b_taxa: taxon_to_group[t] = 'b'
    for t in c_taxa: taxon_to_group[t] = 'c'
    for t in d_taxa: taxon_to_group[t] = 'd'

    return taxon_to_group, a_taxa, b_taxa, c_taxa, d_taxa


# =============================================================================
# Data loading
# =============================================================================

def load_significant_quartets(filename, correction_method, alpha):
    """
    Load significant quartet rows from a .quartet_asym.tsv file.

    Returns a defaultdict mapping
        (node1, node2, quadripartition) -> list of quartet dicts
    containing only rows that pass the significance filter.
    """
    results = defaultdict(list)

    with open(filename, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        fieldnames = reader.fieldnames or []

        # Validate input file
        required = ['taxon_a', 'taxon_b', 'taxon_c', 'taxon_d',
                    'n_ab_cd', 'n_ac_bd', 'n_ad_bc', 'pvalue',
                    'node1', 'node2', 'quadripartition']
        missing = [c for c in required if c not in fieldnames]
        if missing:
            if 'n_quartets' in fieldnames or 'mean_ab_cd' in fieldnames:
                raise ValueError(
                    "'{}' looks like a quadripartition_summary.tsv.\n"
                    "This script requires the detailed .quartet_asym.tsv "
                    "file.".format(filename)
                )
            raise ValueError(
                "Missing required columns: {}".format(', '.join(missing)))

        if correction_method == 'none':
            sig_col = None
        else:
            sig_col = 'sig_{}'.format(correction_method)
            if sig_col not in fieldnames:
                available = [c.replace('sig_', '')
                             for c in fieldnames if c.startswith('sig_')]
                if not available:
                    raise ValueError(
                        "No significance columns found. Expected columns "
                        "like 'sig_BH', 'sig_bonferroni', etc.")
                raise ValueError(
                    "Correction '{}' not found. Available: {}".format(
                        correction_method, ', '.join(available)))

        for row in reader:
            is_sig = (float(row['pvalue']) < alpha
                      if sig_col is None
                      else row[sig_col] == 'yes')
            if is_sig:
                key = (row['node1'], row['node2'], row['quadripartition'])
                results[key].append({
                    'a':       row['taxon_a'],
                    'b':       row['taxon_b'],
                    'c':       row['taxon_c'],
                    'd':       row['taxon_d'],
                    'n_ab_cd': int(row['n_ab_cd']),
                    'n_ac_bd': int(row['n_ac_bd']),
                    'n_ad_bc': int(row['n_ad_bc']),
                    'pvalue':  float(row['pvalue']),
                })

    return results


def load_all_quartet_counts(filename):
    """
    Return a dict mapping (node1, node2, quadripartition) -> total quartet count.
    Used to report the proportion of significant quartets per branch.
    """
    counts = {}
    with open(filename, 'r') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            key = (row['node1'], row['node2'], row['quadripartition'])
            counts[key] = counts.get(key, 0) + 1
    return counts


# =============================================================================
# Direction statistics
# =============================================================================

def compute_direction_stats(quartets):
    """
    Count quartets favouring ac|bd vs ad|bc (ties excluded) and run a
    two-sided binomial test of H0: p(ac|bd favoured) = 0.5.

    Uses log-space arithmetic to avoid integer overflow for large quartet
    counts (naive comb() overflows at ~1,050 evenly-split quartets).

    Returns (n_ac_bd_favoured, n_ad_bc_favoured, binomial_pvalue).
    """
    n_ac = sum(1 for q in quartets if q['n_ac_bd'] > q['n_ad_bc'])
    n_ad = sum(1 for q in quartets if q['n_ad_bc'] > q['n_ac_bd'])
    total = n_ac + n_ad

    if total == 0:
        return n_ac, n_ad, 1.0

    log_half = math.log(0.5)

    def log_binom_pmf(k, n):
        if k < 0 or k > n:
            return -math.inf
        lc = (math.lgamma(n + 1)
              - math.lgamma(k + 1)
              - math.lgamma(n - k + 1))
        return lc + n * log_half   # p = 0.5 → p^k*(1-p)^(n-k) = 0.5^n

    log_obs  = log_binom_pmf(min(n_ac, n_ad), total)
    log_pval = -math.inf
    eps      = 1e-10

    for k in range(total + 1):
        lp = log_binom_pmf(k, total)
        if lp <= log_obs + eps:
            if lp > log_pval:
                log_pval = lp + math.log1p(math.exp(log_pval - lp))
            else:
                log_pval = log_pval + math.log1p(math.exp(lp - log_pval))

    pval = min(math.exp(log_pval) if log_pval > -math.inf else 0.0, 1.0)
    return n_ac, n_ad, pval


# =============================================================================
# Taxon tallying
# =============================================================================

def analyze_branch(sig_quartets):
    """
    Tally taxon appearances and per-taxon ratio statistics among significant
    quartets.

    For each taxon the favored direction is determined by the taxon-level lean
    (whichever of ac_bd_favoured or ad_bc_favoured is larger).  All ratio
    calculations use that taxon-level direction consistently.

    Ratio 1  (r1): n_favored_minority / n_disfavored_minority
      Excluded when n_disfavored_minority == 0; these contribute to r1_excl.

    Ratio 2  (r2): n_favored_minority / n_ab_cd
      Excluded when n_ab_cd == 0; these contribute to r2_excl.

    Count 3  (gt_maj): number of quartets where n_favored_minority > n_ab_cd.

    Returns:
      taxon_counts    Counter  total significant quartets per taxon
      ac_bd_favoured  Counter  quartets where n_ac_bd > n_ad_bc
      ad_bc_favoured  Counter  quartets where n_ad_bc > n_ac_bd
      r1_sums         dict     taxon -> sum of ratio-1 values (for mean)
      r1_counts       dict     taxon -> number of quartets included in r1 mean
      r1_excl         dict     taxon -> number of quartets excluded from r1
      r2_sums         dict     taxon -> sum of ratio-2 values (for mean)
      r2_counts       dict     taxon -> number of quartets included in r2 mean
      r2_excl         dict     taxon -> number of quartets excluded from r2
      gt_maj          dict     taxon -> count where n_favored_minority > n_ab_cd
    """
    taxon_counts   = Counter()
    ac_bd_favoured = Counter()
    ad_bc_favoured = Counter()

    # Accumulate per-quartet data keyed by taxon for ratio calculation later.
    # We store the raw quartet values; ratios are computed after direction is
    # determined from the full pass.
    taxon_quartets = {}   # taxon -> list of quartet dicts it appears in

    for q in sig_quartets:
        for taxon in (q['a'], q['b'], q['c'], q['d']):
            taxon_counts[taxon] += 1
            taxon_quartets.setdefault(taxon, []).append(q)
        if q['n_ac_bd'] > q['n_ad_bc']:
            for taxon in (q['a'], q['b'], q['c'], q['d']):
                ac_bd_favoured[taxon] += 1
        elif q['n_ad_bc'] > q['n_ac_bd']:
            for taxon in (q['a'], q['b'], q['c'], q['d']):
                ad_bc_favoured[taxon] += 1

    # Now compute ratios using taxon-level favored direction
    r1_sums   = {}
    r1_counts = {}
    r1_excl   = {}
    r2_sums   = {}
    r2_counts = {}
    r2_excl   = {}
    gt_maj    = {}

    for taxon, quartets in taxon_quartets.items():
        # Determine taxon-level favored direction
        use_ac = ac_bd_favoured.get(taxon, 0) >= ad_bc_favoured.get(taxon, 0)

        rs1, rc1, re1 = 0.0, 0, 0
        rs2, rc2, re2 = 0.0, 0, 0
        gm = 0

        for q in quartets:
            n_fav  = q['n_ac_bd'] if use_ac else q['n_ad_bc']
            n_dis  = q['n_ad_bc'] if use_ac else q['n_ac_bd']
            n_maj  = q['n_ab_cd']

            # Ratio 1
            if n_dis == 0:
                re1 += 1
            else:
                rs1 += n_fav / n_dis
                rc1 += 1

            # Ratio 2
            if n_maj == 0:
                re2 += 1
            else:
                rs2 += n_fav / n_maj
                rc2 += 1

            # Count 3
            if n_fav > n_maj:
                gm += 1

        r1_sums[taxon]   = rs1
        r1_counts[taxon] = rc1
        r1_excl[taxon]   = re1
        r2_sums[taxon]   = rs2
        r2_counts[taxon] = rc2
        r2_excl[taxon]   = re2
        gt_maj[taxon]    = gm

    return (taxon_counts, ac_bd_favoured, ad_bc_favoured,
            r1_sums, r1_counts, r1_excl,
            r2_sums, r2_counts, r2_excl,
            gt_maj)


# =============================================================================
# Report formatting
# =============================================================================

def _stars(p):
    return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''


def _fmt_ratio(total, ratio_sum, ratio_count, excl, is_r1=True):
    """
    Format a mean ratio with optional exclusion parenthetical.

    If ratio_count == 0 (all excluded): return 'inf(N)' where N = total.
    If excl > 0:                        return 'X.XX(E)' where E = excl count.
    Otherwise:                          return 'X.XX'.

    For ratio 1, all-excluded means every disfavored minority count was zero,
    which represents the strongest possible directional signal (hence inf).
    """
    if ratio_count == 0:
        return 'inf({:d})'.format(total)
    mean = ratio_sum / ratio_count
    base = '{:.2f}'.format(mean)
    return '{:s}({:d})'.format(base, excl) if excl > 0 else base


def format_report(sig_results, title, correction_label, alpha, total_counts=None):
    """
    Format the taxon-involvement text report.

    sig_results : dict mapping (node1, node2, quadripartition) -> list of
                  quartet dicts (only significant ones).
    title       : section heading string.
    """
    lines = []
    lines.append('=' * 78)
    lines.append(title)
    lines.append('=' * 78)
    lines.append('')
    lines.append('Correction method: {}'.format(correction_label))
    lines.append('Alpha level:       {}'.format(alpha))
    lines.append('')
    lines.append('Direction interpretation:')
    lines.append('  Positive bias (+): taxon in quartets where ac|bd > ad|bc')
    lines.append('    (left-side taxa cluster with group_c more than group_d)')
    lines.append('  Negative bias (-): taxon in quartets where ad|bc > ac|bd')
    lines.append('    (left-side taxa cluster with group_d more than group_c)')
    lines.append('')
    lines.append('Column descriptions:')
    lines.append('  ac|bd / ad|bc  Quartets where each minority topology is favored')
    lines.append('  Bias           (ac|bd - ad|bc) / (ac|bd + ad|bc); range [-1, +1]')
    lines.append('  R1             Mean(n_favored_minority / n_disfavored_minority);')
    lines.append('                 direction set by taxon-level lean.')
    lines.append('                 inf(N): all N quartets had disfavored minority = 0')
    lines.append('                 X.XX(N): mean of included quartets; N excluded')
    lines.append('  R2             Mean(n_favored_minority / n_ab_cd);')
    lines.append('                 X.XX(N): N quartets excluded where n_ab_cd = 0')
    lines.append('  >maj           Quartets where n_favored_minority > n_ab_cd')
    lines.append('')
    lines.append('The binomial test P-value tests for a directional bias across')
    lines.append('all significant quartets. Fewer than 6 significant quartets')
    lines.append('gives insufficient power to assess directionality.')
    lines.append('')

    if not sig_results:
        lines.append('No significant asymmetric quartets found.')
        lines.append('=' * 78)
        return '\n'.join(lines)

    sorted_branches = sorted(sig_results.keys(),
                             key=lambda x: -len(sig_results[x]))

    for node1, node2, quadripartition in sorted_branches:
        sig_quartets = sig_results[(node1, node2, quadripartition)]
        n_sig = len(sig_quartets)

        _, a_taxa, b_taxa, c_taxa, d_taxa = parse_quadripartition(quadripartition)
        if not a_taxa:
            continue

        (taxon_counts, ac_fav, ad_fav,
         r1_sums, r1_counts, r1_excl,
         r2_sums, r2_counts, r2_excl,
         gt_maj) = analyze_branch(sig_quartets)
        n_ac, n_ad, dir_pval = compute_direction_stats(sig_quartets)

        lines.append('-' * 78)
        lines.append('Branch: {} <--> {}'.format(node1, node2))
        key = (node1, node2, quadripartition)
        if total_counts and key in total_counts:
            total = total_counts[key]
            lines.append('Significant quartets: {} ({:.1f}% of {:,} total)'.format(
                n_sig, 100.0 * n_sig / total, total))
        else:
            lines.append('Significant quartets: {}'.format(n_sig))
        lines.append('')
        lines.append('Quadripartition: {}'.format(quadripartition))
        lines.append('')

        total_dir = n_ac + n_ad
        if total_dir > 0:
            lines.append(
                'Direction: ac|bd favored in {} ({:.1f}%), '
                'ad|bc in {} ({:.1f}%)'.format(
                    n_ac, 100 * n_ac / total_dir,
                    n_ad, 100 * n_ad / total_dir))
            fmt = '{:.3e}' if dir_pval < 0.001 else '{:.4f}'
            lines.append('Binomial test: P = {} {}'.format(
                fmt.format(dir_pval), _stars(dir_pval)))
            if n_sig < 6:
                lines.append('')
                lines.append(
                    'WARNING: Fewer than 6 significant quartets -- '
                    'insufficient power for directionality test')
        lines.append('')

        # Column header widths:
        # Taxon(20) Count(8) %(8) ac|bd(8) ad|bc(8) Bias(6)
        # R1(12) R2(12) >maj(6)
        col_hdr = ('  {:20s} {:>8s} {:>7s} {:>8s} {:>8s} {:>6s}'
                   ' {:>12s} {:>12s} {:>6s}').format(
            'Taxon', 'Count', '%', 'ac|bd', 'ad|bc', 'Bias',
            'R1', 'R2', '>maj')
        col_div = '  ' + '-' * (len(col_hdr) - 2)

        for grp_label, grp_taxa in [
            ('GROUP A (left side, first)',   a_taxa),
            ('GROUP B (left side, second)',  b_taxa),
            ('GROUP C (right side, first)',  c_taxa),
            ('GROUP D (right side, second)', d_taxa),
        ]:
            present = [t for t in grp_taxa if t in taxon_counts]
            if not present:
                continue
            lines.append('  {}:'.format(grp_label))
            lines.append(col_hdr)
            lines.append(col_div)
            for taxon in sorted(present, key=lambda t: (-taxon_counts[t], t)):
                count = taxon_counts[taxon]
                pct   = 100.0 * count / n_sig
                ac    = ac_fav.get(taxon, 0)
                ad    = ad_fav.get(taxon, 0)
                denom = ac + ad
                bias  = '{:+.2f}'.format((ac - ad) / denom) if denom else ''
                r1 = _fmt_ratio(count,
                                r1_sums.get(taxon, 0.0),
                                r1_counts.get(taxon, 0),
                                r1_excl.get(taxon, 0))
                r2 = _fmt_ratio(count,
                                r2_sums.get(taxon, 0.0),
                                r2_counts.get(taxon, 0),
                                r2_excl.get(taxon, 0),
                                is_r1=False)
                gm = gt_maj.get(taxon, 0)
                lines.append(
                    ('  {:20s} {:>8d} {:>6.1f}% {:>8d} {:>8d} {:>6s}'
                     ' {:>12s} {:>12s} {:>6d}').format(
                        taxon, count, pct, ac, ad, bias, r1, r2, gm))
            lines.append('')

    lines.append('=' * 78)
    return '\n'.join(lines)


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Taxon involvement summary for quartet asymmetry results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Reads the detailed .quartet_asym.tsv from rgc_quartet_asym.py and reports
which taxa are involved in significant asymmetric quartets, with per-taxon
counts and directional bias scores.

Examples:
  python quartet_asym_summary.py -i results.quartet_asym.tsv -o OUTPREFIX
  python quartet_asym_summary.py -i results.quartet_asym.tsv -o OUTPREFIX \\
      --correction bonferroni --alpha 0.01
        '''
    )
    parser.add_argument('--input',  '-i', required=True,
                        help='Input .quartet_asym.tsv file from rgc_quartet_asym.py')
    parser.add_argument('--output', '-o', required=True,
                        help='Output prefix (produces OUTPREFIX.summary.txt)')
    parser.add_argument('--correction', '-c', default='BH',
                        choices=['BH', 'BY', 'bonferroni', 'holm', 'none'],
                        help='Multiple-test correction column to use '
                             '(default: BH)')
    parser.add_argument('--alpha', type=float, default=0.05,
                        help='Significance threshold (default: 0.05)')

    args = parser.parse_args()

    try:
        sig_results   = load_significant_quartets(
            args.input, args.correction, args.alpha)
        total_counts  = load_all_quartet_counts(args.input)
    except (ValueError, OSError) as e:
        print("Error: {}".format(e), file=sys.stderr)
        return 1

    n_sig_q    = sum(len(v) for v in sig_results.values())
    n_branches = len(sig_results)

    if args.correction == 'none':
        title  = ('TAXON INVOLVEMENT IN SIGNIFICANT ASYMMETRIC QUARTETS '
                  '(UNCORRECTED)')
        label  = 'none (uncorrected)'
    else:
        title  = ('TAXON INVOLVEMENT IN SIGNIFICANT ASYMMETRIC QUARTETS '
                  '({} CORRECTED)'.format(args.correction.upper()))
        label  = args.correction

    print("Loaded {:,} significant quartets ({} correction) "
          "across {} branch(es)".format(n_sig_q, label, n_branches))

    report = format_report(sig_results, title, label, args.alpha, total_counts)

    out_path = '{}.summary.txt'.format(args.output)
    with open(out_path, 'w') as f:
        f.write(report)
        f.write('\n')

    print("Summary written to: {}".format(out_path))
    return 0


if __name__ == '__main__':
    sys.exit(main())
