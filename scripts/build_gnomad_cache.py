#!/usr/bin/env python3
"""
build_gnomad_cache.py

Parse each gnomAD per-chromosome VCF, collect all SNPs with AF>=1% into an IntervalTree,
and serialize each tree to cache/chr{N}.tree.pkl with a live progress bar.
"""

import os
import pickle
from cyvcf2 import VCF
from intervaltree import IntervalTree
from tqdm import tqdm

# --- CONFIGURATION ----------------------------------------------------------
GNOMAD_DIR = "/labs/congle/jwhui/gnomad_v4.1/genomes"
CACHE_DIR  = "./cache"
MAF_THRESH = 0.01  # threshold for filtering common SNPs

# ensure cache directory exists
os.makedirs(CACHE_DIR, exist_ok=True)

# gather only the .vcf.bgz filenames
vcf_files = [f for f in os.listdir(GNOMAD_DIR) if f.endswith(".vcf.bgz")]

# loop with progress bar
for fn in tqdm(vcf_files, desc="Building chr trees", unit="chr"):
    chrom = fn.split(".chr")[-1].split(".")[0]
    print(f"\\n[{chrom}] Loading {fn} ...")
    tree = IntervalTree()
    vcf  = VCF(os.path.join(GNOMAD_DIR, fn))

    # iterate all records in this chromosome
    for rec in vcf:
        af = rec.INFO.get("AF", 0.0)
        if af >= MAF_THRESH:
            p0 = rec.POS - 1
            tree.addi(p0, p0+1)

    tree.merge_overlaps()

    out_path = os.path.join(CACHE_DIR, f"chr{chrom}.tree.pkl")
    with open(out_path, "wb") as out:
        pickle.dump(tree, out, protocol=4)

    print(f"[{chrom}] {len(tree)} SNP intervals cached → {out_path}")

print("\\nCache build complete!")