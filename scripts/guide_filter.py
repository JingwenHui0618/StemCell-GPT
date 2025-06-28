#!/usr/bin/env python3
import os, re, pickle, sys
from pyfaidx import Fasta
from intervaltree import IntervalTree

# — config ——————————————————————————————————————————————
REF_FASTA = "/labs/congle/jwhui/hg38.fa"
TREE_DIR  = "/labs/congle/jwhui/scripts/cache/trees"

# load reference (uppercase)
REF = Fasta(REF_FASTA, sequence_always_upper=True)

class PopulationSNPs:
    def __init__(self, tree_dir=TREE_DIR):
        self.tree_dir = tree_dir
        self.trees = {}

    def _load_tree(self, chrom):
        # chrom is already "chr7", "chrX", etc.
        key = chrom[3:]  # drop "chr"
        pkl = os.path.join(self.tree_dir, f"chr{key}.tree.pkl")
        if not os.path.isfile(pkl):
            raise FileNotFoundError(f"Missing SNP tree for {chrom}: {pkl}")
        with open(pkl,"rb") as fh:
            self.trees[key] = pickle.load(fh)

    def has_common_snp(self, chrom, start, end):
        # normalize key
        key = chrom[3:]
        if key not in self.trees:
            self._load_tree(chrom)
        # intervaltree uses 0-based, half-open intervals
        return bool(self.trees[key].overlap(start, end))

class GuideDesigner:
    PAM_RE = re.compile(r'(?=(.{20})GG)')  # lookahead: 20nt then "GG"

    def __init__(self, pop_snp):
        self.pop = pop_snp

    def find_guides(self, chrom, start, end):
        seq = REF[chrom][start:end].seq
        for m in self.PAM_RE.finditer(seq):
            g0 = start + m.start(1)
            g1 = g0 + 20
            yield g0, g1, m.group(1)

    def filter_guides(self, chrom, start, end):
        for g0, g1, prot in self.find_guides(chrom, start, end):
            # tier‐1: no common SNP
            if self.pop.has_common_snp(chrom, g0, g1):
                continue
            yield chrom, g0, g1, prot

def main():
    if len(sys.argv) != 4:
        sys.stderr.write(f"Usage: {sys.argv[0]} <chr|num> <start> <end>\n")
        sys.exit(1)

    raw, start, end = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    # normalize "7" -> "chr7", leave "chr7" alone
    chrom = raw if raw.lower().startswith("chr") else "chr" + raw

    pop = PopulationSNPs()
    gd  = GuideDesigner(pop)
    for c, g0, g1, seq in gd.filter_guides(chrom, start, end):
        print(f"{c}\t{g0}\t{g1}\t{seq}")

if __name__ == "__main__":
    main()