# SDPquartets and RGC asymmetry tests

These programs implement a method for phylogenetic estimation and an introgression test using rare genomic change (RGC) data. The programs are intended to be used with binary RGC data, such as transposable element (TE) insertions.

If you use this software please cite:

Springer MS, Molloy EK, Sloan DB, Simmons MP, Gatesy J. 2020. ILS-Aware analysis of low-homoplasy retroelement insertions: Inference of species trees and introgression using quartets. Journal of Heredity 111:147–168, https://doi.org/10.1093/jhered/esz076

- Springer et al. (2020) describes the SDPquartets method and the quartet asymmetry test used in this software.

(once accepted) Wang N, Kimball RT, Xiao M, Liang B, Braun EL. In review. Rate variation can result in misleading phylogenetic signal in analyses of rare genomic changes.

- Wang et al. (in review) describes these specific implementations of SDPquartets and the automatic selection of quartets induced by the set of quadripartitions in a tree used for the asymmetry test.

Swofford DL (2003) PAUP\*. Phylogenetic Analysis Using Parsimony (\*and Other Methods). https://paup.phylosolutions.com

## Requirements

- Python 3.8 or later (standard library only; no third-party packages required).
- PAUP* (Swofford 2003, https://paup.phylosolutions.com) is required by `sdp_quartets.py` for the MRP tree search and bootstrap. The asymmetry-test programs (`rgc_quartet_asym.py`, `quartet_asym_summary.py`) have no external dependencies.

## Input data

All programs read a binary character matrix in relaxed PHYLIP or NEXUS format (the format is auto-detected from the first line). Character states are `0` and `1`; missing data may be coded as `-` or `?`. Only "two against two" site patterns are phylogenetically informative for a quartet; constant sites, singletons (three-against-one), and patterns with missing data are ignored. The asymmetry test additionally requires a species tree in Newick format.

## Species tree estimation with SDPquartets

`sdp_quartets.py` is a python implementation of the SDPquartets (Split Decomposition with Parsimony for Quartets) method proposed by Springer et al. (2020). SDPquartets is a "species tree" method for binary RGCs.

The SDPquartets method:

1. Enumerates all C(n,4) quartets from a binary character matrix
2. Determines the most parsimonious (MP) quartet tree(s) by pattern counting
3. Generates an MRP (Matrix Representation using Parsimony) matrix for the MP quartets
4. Infers the species tree by MP tree for the quartet MRP matrix

The original Perl implementation of SDPquartets is available from https://github.com/dbsloan/SDPquartets

Like the earlier implementation it uses PAUP* (Swofford 2003, https://paup.phylosolutions.com) for the MRP analysis that is used to combine trees, although other parsimony programs could be used for that step. 

The most parsimonious quartets given the binary RGC data are identified within the python code and used to generate the MRP file that PAUP* uses to combine the quartets into a single file.

```
usage: sdp_quartets.py [-h] --input INPUT --output OUTPUT [--format {phylip,nexus}] [--springer-compat] [--paup] [--paup-path PAUP_PATH]
                       [--bootstrap] [--reps REPS] [--seed SEED] [--boot-no-opt] [--lazy] [--clean]
SDPquartets: Split Decomposition with Parsimony for Quartets
options:
  -h, --help            show this help message and exit
  --input INPUT, -i INPUT
                        Input character matrix (PHYLIP or NEXUS format)
  --output OUTPUT, -o OUTPUT
                        Base name for output files
  --format {phylip,nexus}, -f {phylip,nexus}
                        Output MRP matrix format (default: nexus)
  --springer-compat     Include 3-way ties with 6x/3x/2x weighting for compatibility with Springer et al. 2020
  --paup                Run PAUP* tree search on MRP matrix
  --paup-path PAUP_PATH
                        Path to PAUP* executable (default: paup)
  --bootstrap, -b       Perform bootstrap analysis
  --reps REPS           Number of bootstrap replicates (default: 100)
  --seed SEED           Random number seed for bootstrap (default: from system clock)
  --boot-no-opt         Skip optimal tree search when bootstrapping (only do bootstrap)
  --lazy                Use faster but less thorough PAUP* search settings
  --clean               Remove bootstrap directory and intermediate files after analysis
```

### Output files

- `<output>.quartet_log.tsv` — one row per quartet, giving the resolution status (1 = resolved, 2 = two-way tie, 3 = three-way tie) and the parsimony score of each of the three topologies.
- `<output>.quartets.tre` — the MP quartet trees in Newick format.
- `<output>.mrp.nex` or `<output>.mrp.phy` — the MRP matrix (NEXUS by default; PHYLIP with `-f phylip`).
- `<output>.SDPtrees.tre` and `<output>.SDPconstree.tre` — the MP trees and strict consensus from the PAUP* search (with `--paup`).
- `<output>.bootstrap_trees.tre` plus a `<output>/` directory of per-replicate files (with `--bootstrap`).

### MRP encoding and bootstrap weighting

By default, each resolved quartet contributes two identical MRP characters and each two-way tie contributes one character per tied topology (a 2:1 weighting of resolved versus tied quartets); three-way ties are omitted. Adding `--springer-compat` instead applies the 6×/3×/2× weighting of Springer et al. (2020), which retains three-way ties.

With `--bootstrap`, characters in the original binary RGC data matrix are resampled with replacement for each replicate, then quartet analyses are conducted, an MRP dataset is generated, and the PAUP* search is conducted using the MRP data. Trees from a replicate that recovers *n* equally parsimonious trees are each written with weight 1/*n* in `<output>.bootstrap_trees.tre`. 

To make the bootstrap consensus respect those weights, run `contree` with `usetreewts=yes` in PAUP\* (otherwise replicates with tied trees are over-counted). If PAUP\* is not used for the search, the MRP matrix can be handed to any other parsimony program.

## Introgression testing with the quartet asymmetry test

`rgc_quartet_asym.py` implements the quartet asymmetry test of Springer et al. (2020). Given a binary RGC matrix and a species tree, it asks, for each internal branch of the tree, whether the two minority quartet resolutions occur with equal frequency. Under incomplete lineage sorting (ILS) without introgression the two minority resolutions are expected to be equally frequent, so a significant asymmetry is evidence of gene flow.

The method:

1. Reads the species tree and identifies the quadripartition induced by each internal branch (the four groups of taxa surrounding that branch).
2. For each quadripartition, automatically selects the set of valid quartets — one taxon drawn from each of the four groups — so that the species-tree-concordant resolution is `ab|cd` by construction. This is the automatic selection of quartets induced by the quadripartitions described by Wang et al. (in review), and it prevents quartets from being mis-attributed to the wrong branch.
3. Counts the binary site patterns supporting each of the three quartet resolutions (`ab|cd`, `ac|bd`, `ad|bc`).
4. Tests the two minority counts (`ac|bd` versus `ad|bc`) against a 1:1 null with a two-sided exact binomial test.
5. Applies multiple test corrections (Benjamini–Hochberg, Benjamini–Yekutieli, Bonferroni, Holm–Bonferroni) independently within each quadripartition, since each branch is a separate family of hypotheses.

**Note that this program automatically adds numerial indices to the input species tree.** These indices are used to identify branches associated with each quadripartition. If you are using your own branch indices you will have to reconcile the indices used by this program with those you are using. **A newick tree with the indices is part of the output** (see below).

```
usage: rgc_quartet_asym.py [-h] --input INPUT --tree TREE --output OUTPUT
                           [--alpha ALPHA]
                           [--correction {BH,BY,bonferroni,holm} [{BH,BY,bonferroni,holm} ...]]
Quartet Asymmetry Test for Binary Characters (RGC/Indel data)
options:
  -h, --help            show this help message and exit
  --input INPUT, -i INPUT
                        Input character matrix (PHYLIP or NEXUS format)
  --tree TREE, -t TREE  Species tree in Newick format
  --output OUTPUT, -o OUTPUT
                        Base name for output files
  --alpha ALPHA         Significance level (default: 0.05)
  --correction {BH,BY,bonferroni,holm} [{BH,BY,bonferroni,holm} ...]
                        Multiple test correction methods (default: BH BY bonferroni holm)
```

### Output files

- `<output>.quartet_asym.tsv` — detailed per-quartet results. One row per valid quartet, giving the four taxa, the pattern counts (`n_ab_cd`, `n_ac_bd`, `n_ad_bc`), whether `ab|cd` is the majority resolution, the uncorrected binomial p-value, the adjusted p-value and significance flag for each requested correction, and the branch (`node1`, `node2`) and quadripartition the quartet belongs to.
- `<output>.quadripartition_summary.tsv` — one row per internal branch, giving mean pattern counts and the number of significant quartets under each correction.
- `<output>.labeled_tree.tre` — the input species tree with internal nodes labeled by the numeric IDs used in the `node1`/`node2` columns of the other output files.

### Summarizing taxon involvement in asymmetric quartets

`quartet_asym_summary.py` reads the detailed `.quartet_asym.tsv` file and, for each branch with at least one significant quartet, reports the direction of the asymmetry and which taxa drive it. For each branch it gives the overall direction (`ac|bd`- versus `ad|bc`-favored) with a binomial test for directional bias, and, per taxon, how often the taxon appears in significant quartets, a directional bias score, and ratio statistics comparing the favored minority resolution to the disfavored minority and to the majority resolution. This is useful for identifying the lineages involved in a putative introgression event. The output is a single text report, `OUTPREFIX.summary.txt`.

By default the summary uses the Benjamini–Hochberg significance column; use `-c` to select a different correction (matching one used when running `rgc_quartet_asym.py`), or `-c none` to threshold on uncorrected p-values.

```
usage: quartet_asym_summary.py [-h] --input INPUT --output OUTPUT
                               [--correction {BH,BY,bonferroni,holm,none}]
                               [--alpha ALPHA]
Taxon involvement summary for quartet asymmetry results
options:
  -h, --help            show this help message and exit
  --input INPUT, -i INPUT
                        Input .quartet_asym.tsv file from rgc_quartet_asym.py
  --output OUTPUT, -o OUTPUT
                        Output prefix (produces OUTPREFIX.summary.txt)
  --correction {BH,BY,bonferroni,holm,none}, -c {BH,BY,bonferroni,holm,none}
                        Multiple-test correction column to use (default: BH)
  --alpha ALPHA         Significance threshold (default: 0.05)
```

## Example workflow

```
# 1. Estimate a species tree from the RGC matrix
python sdp_quartets.py -i data.phy -o sdp_result --paup

# 2. Test for introgression along the branches of a species tree
python rgc_quartet_asym.py -i data.phy -t species.tre -o asym_result

# 3. Summarize which taxa drive the significant asymmetries
python quartet_asym_summary.py -i asym_result.quartet_asym.tsv -o asym_result
```

## Test datasets

The `SDPasym_TEST_DATA/` directory bundles input data and example output for both programs, using two of the empirical data sets from Springer et al. (2020): **Balaenopteroidea** (five rorquals and the gray whale plus a right-whale outgroup; 7 taxa, 17,225 retroelement characters) and **Palaeognathae** (12 ratites and tinamous plus a chicken outgroup; 13 taxa, 4,301 retroelement characters).

Although the asymmetry tests can use any estimate of the species tree as input, including the SDPquartets trees, the asymmetry test output in these folder use ASTRAL_BP estimates of the species trees. Those trees were distributed in the supporting information of Springer et al. (2020) and they are included in the `00_Data_files_for_test` folder.

### Folder contents

- `00_Data_files_for_test/` — the binary RGC matrices (`*.nex`) and the ASTRAL_BP species trees (`*_ASTRAL_BP_species.tre`).
- `02_Balaenopteroidea_SDP_test/` and `03_Palaeognathae_SDP_test/` — output of `sdp_quartets.py` (PAUP* search plus 100 bootstrap replicates): the run log (`*_RUN_TEST.txt`), the MRP matrix, the MP and consensus trees, and the weighted bootstrap treefile. The bootstrap consensus in each log was computed with `contree all / majrule=yes le50=yes usetreewts=yes;` so that the 1/*n* tree weights are applied.
- `04_Asymmetry_result_files/` — output of `rgc_quartet_asym.py` and `quartet_asym_summary.py` for both data sets: the per-quartet table, the quadripartition summary, the labeled tree, a run log, and the taxon-involvement summary.

### Expected results

These data sets reproduce the main qualitative findings of Springer et al. (2020): the quartet asymmetry test detects significant minority-quartet asymmetry (evidence of introgression) in Balaenopteroidea, but not in Palaeognathae after multiple-test correction.

### Reproducing the analyses

Run from inside `00_Data_files_for_test/` (exact bootstrap values require the random seed recorded in each SDP run log):

```
# Balaenopteroidea
python sdp_quartets.py -i Balaenopteroidea.nex -o Balaenopteroidea_SDP --paup --bootstrap --reps 100
python rgc_quartet_asym.py -i Balaenopteroidea.nex -t Balaenopteroidea_ASTRAL_BP_species.tre -o QUARTETasym_Balaenopteroidea
python quartet_asym_summary.py -i QUARTETasym_Balaenopteroidea.quartet_asym.tsv -o QUARTETasym_Balaenopteroidea

# Palaeognathae
python sdp_quartets.py -i Palaeognathae.nex -o Palaeognathae_SDP --paup --bootstrap --reps 100
python rgc_quartet_asym.py -i Palaeognathae.nex -t Palaeognathae_ASTRAL_BP_species.tre -o QUARTETasym_Palaeognathae
python quartet_asym_summary.py -i QUARTETasym_Palaeognathae.quartet_asym.tsv -o QUARTETasym_Palaeognathae
```

The bootstrap consensus for each SDP run is then computed in PAUP* on `<output>.bootstrap_trees.tre` with `contree all / majrule=yes le50=yes usetreewts=yes;`.

### Taxon abbreviations

Both data sets label taxa with short codes; the outgroup is listed last in each table.

**Balaenopteroidea**

| Code | Common name | Species |
|------|-------------|---------|
| MW | minke whale | *Balaenoptera acutorostrata* |
| GW | gray whale | *Eschrichtius robustus* |
| FW | fin whale | *Balaenoptera physalus* |
| HW | humpback whale | *Megaptera novaeangliae* |
| SW | sei whale | *Balaenoptera borealis* |
| BW | blue whale | *Balaenoptera musculus* |
| NA | North Atlantic right whale (outgroup) | *Eubalaena glacialis* |

**Palaeognathae**

| Code | Common name | Species |
|------|-------------|---------|
| aptHaa | great spotted kiwi | *Apteryx haastii* |
| aptOwe | little spotted kiwi | *Apteryx owenii* |
| aptRow | Okarito brown kiwi | *Apteryx rowi* |
| casCas | southern cassowary | *Casuarius casuarius* |
| droNov | emu | *Dromaius novaehollandiae* |
| rheAme | greater rhea | *Rhea americana* |
| rhePen | lesser (Darwin's) rhea | *Rhea pennata* |
| cryCin | thicket tinamou | *Crypturellus cinnamomeus* |
| tinGut | white-throated tinamou | *Tinamus guttatus* |
| eudEle | elegant crested tinamou | *Eudromia elegans* |
| notPer | Chilean tinamou | *Nothoprocta perdicaria* |
| strCam | ostrich | *Struthio camelus* |
| Gallus | chicken (outgroup) | *Gallus gallus* |

