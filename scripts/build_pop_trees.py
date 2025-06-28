#!/usr/bin/env python3
import os, pickle
from intervaltree import IntervalTree

# where your .pos.txt files live:
POS_DIR = "cache/snps_pos"
# where we’ll write the pickles:
OUT_DIR = "cache/trees"
os.makedirs(OUT_DIR, exist_ok=True)

def build_chr_tree(chrom):
    txt = os.path.join(POS_DIR, f"chr{chrom}.pos.txt")
    tree = IntervalTree()
    total = 0
    with open(txt) as fh:
        for line in fh:
            pos = int(line.strip())
            # interval [pos-1,pos)
            tree.addi(pos-1, pos)
            total += 1
            if total % 100_000 == 0:
                print(f"  chr{chrom}: added {total} SNPs…")
    tree.merge_overlaps()
    out = os.path.join(OUT_DIR, f"chr{chrom}.tree.pkl")
    with open(out, "wb") as fo:
        pickle.dump(tree, fo)
    print(f"→ chr{chrom} done: {len(tree)} intervals stored in {out}")

if __name__ == "__main__":
    chroms = [str(i) for i in range(1,23)] + ["X","Y"]
    for c in chroms:
        print(f"Building tree for chr{c} …")
        build_chr_tree(c)