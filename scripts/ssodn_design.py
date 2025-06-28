import argparse
from Bio.Seq import reverse_complement
from Bio.Data import CodonTable

# Standard codon table for synonymous mutation lookup
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]


def introduce_silent_mutations(seq, pam_start, edit_rel):
    """
    Introduce silent mutations codon-by-codon from the PAM region until the edit site.
    Returns updated sequence and list of mutation annotations.
    """
    annotations = []
    # Align to codon boundary
    pos = pam_start - (pam_start % 3)
    while pos <= edit_rel and pos + 3 <= len(seq):
        codon = seq[pos:pos+3]
        aa = standard_table.forward_table.get(codon)
        if aa and not (pos <= edit_rel < pos+3):
            synonyms = [c for c,a in standard_table.forward_table.items() if a==aa and c!=codon]
            if synonyms:
                new_codon = synonyms[0]
                seq = seq[:pos] + new_codon + seq[pos+3:]
                annotations.append({'pos': pos, 'from': codon, 'to': new_codon})
        pos += 3
    return seq, annotations


def design_ssodn(hg, chrom, sgRNA,
                  edit_type, edit_coord, edit_to,
                  arm_length=45,
                  start=None, end=None):
    """
    Designs an ssODN donor sequence.
    If start/end provided, restricts guide search to that region; otherwise scans full chromosome.
    Validates that the edit coordinate falls within the homology arms around the sgRNA cut site.
    Supports substitution, insertion, deletion, plus silent PAM-site mutations.
    """
    strand = '+'
    # fetch region or full chromosome
    if start is not None and end is not None:
        seq_region = hg[chrom][start:end].seq.upper()
    else:
        seq_region = hg[chrom][:].seq.upper()
        start = 0

    # locate guide on region
    idx = seq_region.find(sgRNA)
    if idx < 0:
        rc = reverse_complement(sgRNA)
        idx = seq_region.find(rc)
        if idx < 0:
            raise RuntimeError("sgRNA not found in target sequence")
        strand = '-'
        sgRNA = rc

    # genomic cut position (3bp upstream of PAM)
    cut = start + idx + len(sgRNA) - 3

    # validate edit is within homology arms
    if not (cut - arm_length <= edit_coord < cut + arm_length):
        raise ValueError(
            "edit_coord (%d) must be within ±%dbp of cut site at %d" % (edit_coord, arm_length, cut))

    # build homology arms
    left_start = max(0, cut - arm_length)
    left_arm = hg[chrom][left_start:cut].seq.upper()
    right_arm = hg[chrom][cut:cut + arm_length].seq.upper()
    donor = left_arm + right_arm

    ann = {'strand': strand,
           'guide_index': (start + idx, start + idx + len(sgRNA)),
           'cut_site': cut,
           'arms': {'left': (left_start, cut), 'right': (cut, cut + arm_length)}}

    # relative position for edit in donor
    if edit_coord < cut:
        rel = edit_coord - left_start
    else:
        rel = len(left_arm) + (edit_coord - cut)

    # apply edit
    if edit_type == 'substitution':
        donor = donor[:rel] + edit_to + donor[rel+1:]
        ann['edit'] = {'type': 'sub', 'coord': edit_coord, 'new': edit_to}
    elif edit_type == 'insertion':
        donor = donor[:rel] + edit_to + donor[rel:]
        ann['edit'] = {'type': 'ins', 'coord': edit_coord, 'seq': edit_to}
    elif edit_type == 'deletion':
        length = len(edit_to)
        donor = donor[:rel] + donor[rel+length:]
        ann['edit'] = {'type': 'del', 'coord': edit_coord, 'length': length}
    else:
        raise ValueError("Unsupported edit type: %s" % edit_type)

    # silent mutations around PAM until edit site
    pam_rel = len(left_arm) + len(sgRNA) - 3
    donor, silent_ann = introduce_silent_mutations(donor, pam_rel, rel)
    if silent_ann:
        ann['silent_mutations'] = silent_ann

    return {'donor': donor, 'annotations': ann}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ssODN design module")
    parser.add_argument('--chrom', required=True)
    parser.add_argument('--sgRNA', required=True)
    parser.add_argument('--edit_type', choices=['substitution','insertion','deletion'], required=True)
    parser.add_argument('--edit_coord', type=int, required=True)
    parser.add_argument('--edit_to', required=True)
    parser.add_argument('--arm_length', type=int, default=45)
    parser.add_argument('--start', type=int, help='start pos of search region')
    parser.add_argument('--end', type=int, help='end pos of search region')
    args = parser.parse_args()

    import pyfaidx
    hg = pyfaidx.Fasta("/labs/congle/jwhui/hg38.fa")

    result = design_ssodn(
        hg, args.chrom, args.sgRNA,
        args.edit_type, args.edit_coord,
        args.edit_to, args.arm_length,
        args.start, args.end
    )
    print(result)