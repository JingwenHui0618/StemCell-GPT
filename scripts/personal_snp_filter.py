#!/usr/bin/env python3
"""
personal_snp_filter.py

A script to generate CRISPR SpCas9 guide RNAs (20-nt protospacers + NGG PAM)
within a given region, filtering based on overlap with personal SNPs
from a bgzip–compressed, tabix–indexed VCF.
"""
import argparse
import bisect
from pyfaidx import Fasta
from cyvcf2 import VCF

def load_personal_variants(vcf_path):
    variants = {}
    vcf = VCF(vcf_path)
    for var in vcf:
        chrom, pos = var.CHROM, var.POS
        gt0, gt1 = var.genotypes[0][0], var.genotypes[0][1]
        zyg = 'hom' if gt0 == gt1 else 'het'
        for alt in var.ALT:
            variants.setdefault(chrom, []).append((pos, var.REF, alt, zyg))
    for chrom in variants:
        variants[chrom].sort(key=lambda x: x[0])
    return variants

def has_overlap(variants, chrom, start, end):
    hits = []
    if chrom in variants:
        positions = [v[0] for v in variants[chrom]]
        i = bisect.bisect_left(positions, start)
        while i < len(positions) and positions[i] < end:
            hits.append(variants[chrom][i])
            i += 1
    return hits

def find_spcas9_guides(hg, chrom, start, end):
    seq = hg[chrom][start:end].seq.upper()
    guides = []
    for i in range(len(seq)-2):
        if seq[i+1:i+3] == 'GG':
            ps = start + i - 20
            if ps >= start:
                g = hg[chrom][ps:ps+20].seq.upper()
                if len(g)==20 and 'N' not in g:
                    guides.append((ps, '+', g))
        if seq[i:i+2] == 'CC':
            pam_end = start + i + 3
            if pam_end + 20 <= end:
                sub = hg[chrom][pam_end:pam_end+20]
                g = str(sub.reverse.complement).upper()
                if len(g)==20 and 'N' not in g:
                    guides.append((pam_end, '-', g))
    return guides

def main():
    parser = argparse.ArgumentParser(
        description="Personal SNP-aware SpCas9 guide design"
    )
    parser.add_argument('--vcf',        required=True,
                        help='bgzip-compressed, indexed VCF')
    parser.add_argument('--region',     required=True,
                        help='chr:start-end, e.g. chr7:5501900-5503000')
    parser.add_argument('--fasta',
                        default='/labs/congle/jwhui/hg38.fa',
                        help='reference FASTA path')
    parser.add_argument('--out',
                        default='personal_guides.tsv',
                        help='output TSV filename')
    parser.add_argument('--overlap-only', action='store_true',
                        help='only report guides that overlap personal SNPs')
    args = parser.parse_args()

    chrom, coords = args.region.split(':')
    start, end = map(int, coords.split('-'))

    print('Loading reference…')
    hg = Fasta(args.fasta)
    print('Loading personal variants…')
    variants = load_personal_variants(args.vcf)
    print('Discovering guides…')
    guides = find_spcas9_guides(hg, chrom, start, end)
    print('Filtering and writing output…')
    with open(args.out, 'w') as w:
        w.write('start\tstrand\tguide\toverlaps\n')
        for ps, strand, seq in guides:
            hits = has_overlap(variants, chrom, ps, ps+20)
            if args.overlap_only:
                if hits:
                    strs = [f"{p}:{r}>{a}({z})" for p,r,a,z in hits]
                    w.write(f"{ps}\t{strand}\t{seq}\t{';'.join(strs)}\n")
            else:
                if not hits:
                    w.write(f"{ps}\t{strand}\t{seq}\tNA\n")
    print(f"Done. Results written to {args.out}")

if __name__ == '__main__':
    main()