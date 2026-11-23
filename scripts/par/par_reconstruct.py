#!/usr/bin/env python3
"""Map PostprocessGermlineCNVCalls output from pseudo-contigs back onto chrX.

Applies to genotyped_segments.vcf.gz, genotyped_intervals.vcf.gz and
denoised_copy_ratios.tsv of one sample (each is optional; give input and output together).

For every file:
  * chrPAR1 / chrPAR2 records are renamed to chrX (coordinates are already original);
  * the VCF '##contig' lines / TSV '@SQ' lines are replaced by the ORIGINAL dictionary
    (the pseudo dictionary minus chrPAR1/chrPAR2, or --original-dict);
  * records are sorted in original dictionary order (PAR1 first on chrX, PAR2 last).

VCFs additionally get INFO tags
  PAR=PAR1|PAR2           PAR records only
  REFCN=<int>             2 on PAR records; on chrX records the sample's chrX ploidy from
                          --contig-ploidy (tag omitted when no ploidy file is given)
  PAR_BOUNDARY_ADJACENT   flag, segments VCF only, on NON-REFERENCE segments touching a
                          PAR / X-specific boundary (see below)

PAR_BOUNDARY_ADJACENT definition. Boundary positions are taken from the sample's intervals
VCF (--intervals-vcf, pseudo-contig version), i.e. from the intervals that were actually
modelled:
  PAR1 side   non-ref chrPAR1 segment with END >= end of the last chrPAR1 interval
  X side      non-ref chrX segment with POS <= start of the first chrX interval
              (the first X-specific interval after PAR1)
  PAR2 side   non-ref chrPAR2 segment with POS <= start of the first chrPAR2 interval
  X side      non-ref chrX segment with END >= end of the last chrX interval
              (the last X-specific interval before PAR2)
--boundary-tolerance (bp, default 0) loosens each comparison. Without --intervals-vcf the
nominal PAR coordinates are used instead (then a tolerance of about one bin is advisable).
Segments are never merged.

Non-reference means ALT != '.' and a genotype that is not all-0 (same logic as the
finalize_genotyped_segments rule).

The VCF output is plain text; sort/compress/index it afterwards with
'bcftools sort -O z' and 'bcftools index -t' (the Snakemake rule does this).
"""

import argparse
import sys

import par_common as pc

INFO_HEADER_LINES = {
    "PAR": '##INFO=<ID=PAR,Number=1,Type=String,Description="Record lies in the pseudo-autosomal '
           'region PAR1 or PAR2 (modelled as diploid pseudo-contig, mapped back to chrX)">',
    "REFCN": '##INFO=<ID=REFCN,Number=1,Type=Integer,Description="Reference copy number used as '
             'baseline for this record (2 in PAR, chrX ploidy call elsewhere on chrX)">',
    "PAR_BOUNDARY_ADJACENT": '##INFO=<ID=PAR_BOUNDARY_ADJACENT,Number=0,Type=Flag,Description='
                             '"Non-reference segment touches a PAR / X-specific boundary; the '
                             'event may continue on the other side, review as one event">',
}


def original_dictionary(pseudo_dict_path, original_dict_path=None):
    """Return original @SQ lines."""
    if original_dict_path:
        return pc.sq_lines_of(pc.read_header(original_dict_path))
    return pc.original_sq_lines(pc.sq_lines_of(pc.read_header(pseudo_dict_path)))


# --------------------------------------------------------------------------------------
# Boundaries
# --------------------------------------------------------------------------------------

def nominal_boundaries():
    par1_end = pc.PAR_COORDS["chrX"]["PAR1"][1]
    par2_start = pc.PAR_COORDS["chrX"]["PAR2"][0]
    return {"par1_last_end": par1_end, "x_first_start": par1_end + 1,
            "par2_first_start": par2_start, "x_last_end": par2_start - 1}


def scan_boundaries(intervals_vcf):
    """Boundary positions from the modelled intervals (pseudo-contig intervals VCF)."""
    par1_end = pc.PAR_COORDS["chrX"]["PAR1"][1]
    par2_start = pc.PAR_COORDS["chrX"]["PAR2"][0]
    b = {"par1_last_end": None, "x_first_start": None,
         "par2_first_start": None, "x_last_end": None}
    with pc.open_text(intervals_vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            contig, pos = f[0], int(f[1])
            end = pc.info_end(f[7], pos)
            if contig == "chrPAR1":
                if b["par1_last_end"] is None or end > b["par1_last_end"]:
                    b["par1_last_end"] = end
            elif contig == "chrPAR2":
                if b["par2_first_start"] is None or pos < b["par2_first_start"]:
                    b["par2_first_start"] = pos
            elif contig == "chrX":
                if pos > par1_end and (b["x_first_start"] is None or pos < b["x_first_start"]):
                    b["x_first_start"] = pos
                if end < par2_start and (b["x_last_end"] is None or end > b["x_last_end"]):
                    b["x_last_end"] = end
    return b


def touches_boundary(contig, pos, end, b, tol):
    if contig == "chrPAR1":
        return b["par1_last_end"] is not None and end >= b["par1_last_end"] - tol
    if contig == "chrPAR2":
        return b["par2_first_start"] is not None and pos <= b["par2_first_start"] + tol
    if contig == "chrX":
        if b["x_first_start"] is not None and pos <= b["x_first_start"] + tol:
            return True
        return b["x_last_end"] is not None and end >= b["x_last_end"] - tol
    return False


# --------------------------------------------------------------------------------------
# VCF
# --------------------------------------------------------------------------------------

def reconstruct_vcf(in_path, out_path, orig_sq, boundaries=None, chrx_ploidy=None,
                    segments=False, tolerance=0):
    """Rewrite one pseudo-contig VCF onto chrX. Returns number of records written."""
    lengths = pc.contig_lengths(orig_sq)
    order = pc.contig_order(orig_sq)
    meta, column_line, records = [], None, []
    with pc.open_text(in_path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith("##"):
                meta.append(line)
            elif line.startswith("#"):
                column_line = line
            elif line:
                records.append(line)
    if column_line is None:
        raise ValueError("%s: no #CHROM header line" % in_path)

    # header: original contigs at the position of the first ##contig line; new INFO lines last
    contig_lines = ["##contig=<ID=%s,length=%d>" % (name, ln) for name, ln in lengths.items()]
    new_meta, placed = [], False
    for line in meta:
        if line.startswith("##contig="):
            if not placed:
                new_meta.extend(contig_lines)
                placed = True
        else:
            new_meta.append(line)
    if not placed:
        new_meta.extend(contig_lines)
    wanted = ["PAR", "REFCN"] + (["PAR_BOUNDARY_ADJACENT"] if segments else [])
    present = [m for m in new_meta if m.startswith("##INFO=")]
    for tag in wanted:
        if not any(m.startswith("##INFO=<ID=%s," % tag) for m in present):
            new_meta.append(INFO_HEADER_LINES[tag])

    out_records = []
    for line in records:
        f = line.split("\t")
        contig, pos = f[0], int(f[1])
        end = pc.info_end(f[7], pos)
        tags = []
        if contig in pc.REGION_FOR_PSEUDO:
            tags.append("PAR=" + pc.REGION_FOR_PSEUDO[contig])
            tags.append("REFCN=2")
        elif contig == "chrX" and chrx_ploidy is not None:
            tags.append("REFCN=%d" % chrx_ploidy)
        if (segments and boundaries is not None and pc.is_non_reference(f)
                and touches_boundary(contig, pos, end, boundaries, tolerance)):
            tags.append("PAR_BOUNDARY_ADJACENT")
        if tags:
            f[7] = ";".join(tags) if f[7] == "." else f[7] + ";" + ";".join(tags)
        if contig in pc.REGION_FOR_PSEUDO:
            f[0] = contig = "chrX"
        if contig not in order:
            raise ValueError("%s: contig %r not in the original dictionary" % (in_path, contig))
        out_records.append((order[contig], pos, "\t".join(f)))
    out_records.sort(key=lambda r: (r[0], r[1]))

    with pc.open_text(out_path, "w") as out:
        for line in new_meta:
            out.write(line + "\n")
        out.write(column_line + "\n")
        for _, _, line in out_records:
            out.write(line + "\n")
    return len(out_records)


# --------------------------------------------------------------------------------------
# denoised copy ratios TSV
# --------------------------------------------------------------------------------------

def reconstruct_tsv(in_path, out_path, orig_sq):
    header, column_header, rows = pc.load_table(in_path)
    order = pc.contig_order(orig_sq)
    new_header = pc.replace_sq_block(header, orig_sq)
    out_rows = []
    for line in rows:
        contig = line.split("\t", 1)[0]
        if contig in pc.REGION_FOR_PSEUDO:
            line = "chrX" + line[len(contig):]
            contig = "chrX"
        if contig not in order:
            raise ValueError("%s: contig %r not in the original dictionary" % (in_path, contig))
        out_rows.append((order[contig], int(line.split("\t", 2)[1]), line))
    out_rows.sort(key=lambda r: (r[0], r[1]))
    pc.write_table(out_path, new_header, column_header, [r[2] for r in out_rows])
    return len(out_rows)


# --------------------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pseudo-dict", required=True,
                    help="pseudo .dict from par_make_dict.py (original dictionary is derived from it)")
    ap.add_argument("--original-dict", default=None,
                    help="optional explicit original dictionary instead of deriving it")
    ap.add_argument("--contig-ploidy", default=None, help="sample's contig_ploidy.tsv (for REFCN on chrX)")
    ap.add_argument("--segments-vcf")
    ap.add_argument("--out-segments-vcf")
    ap.add_argument("--intervals-vcf")
    ap.add_argument("--out-intervals-vcf")
    ap.add_argument("--denoised-copy-ratios")
    ap.add_argument("--out-denoised-copy-ratios")
    ap.add_argument("--boundary-tolerance", type=int, default=0,
                    help="bp tolerance for PAR_BOUNDARY_ADJACENT (default 0)")
    args = ap.parse_args(argv)

    pairs = [("--segments-vcf", args.segments_vcf, args.out_segments_vcf),
             ("--intervals-vcf", args.intervals_vcf, args.out_intervals_vcf),
             ("--denoised-copy-ratios", args.denoised_copy_ratios, args.out_denoised_copy_ratios)]
    for name, src, dst in pairs:
        if bool(src) != bool(dst):
            ap.error("%s and its --out- counterpart must be given together" % name)
    if not any(src for _, src, _ in pairs):
        ap.error("nothing to do: give at least one input/output pair")

    orig_sq = original_dictionary(args.pseudo_dict, args.original_dict)
    chrx_ploidy = None
    if args.contig_ploidy:
        _, ploidy = pc.read_contig_ploidy(args.contig_ploidy)
        if "chrX" in ploidy:
            chrx_ploidy = ploidy["chrX"][0]
        else:
            sys.stderr.write("par_reconstruct: chrX not in contig ploidy calls; REFCN omitted on chrX\n")

    if args.intervals_vcf:
        boundaries = scan_boundaries(args.intervals_vcf)
    else:
        boundaries = nominal_boundaries()

    if args.segments_vcf:
        n = reconstruct_vcf(args.segments_vcf, args.out_segments_vcf, orig_sq, boundaries,
                            chrx_ploidy, segments=True, tolerance=args.boundary_tolerance)
        sys.stderr.write("par_reconstruct: %d segment records\n" % n)
    if args.intervals_vcf:
        n = reconstruct_vcf(args.intervals_vcf, args.out_intervals_vcf, orig_sq, None,
                            chrx_ploidy, segments=False)
        sys.stderr.write("par_reconstruct: %d interval records\n" % n)
    if args.denoised_copy_ratios:
        n = reconstruct_tsv(args.denoised_copy_ratios, args.out_denoised_copy_ratios, orig_sq)
        sys.stderr.write("par_reconstruct: %d denoised copy ratio rows\n" % n)
    return 0


if __name__ == "__main__":
    sys.exit(main())
