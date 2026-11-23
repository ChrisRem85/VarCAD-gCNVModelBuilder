"""Tests for scripts/par on small synthetic files (no GATK needed)."""

import gzip
import os
import subprocess
import sys

import pytest

import par_common as pc
import par_make_dict
import par_ploidy_qc
import par_reconstruct
import par_relabel

SCRIPTS = os.path.join(os.path.dirname(__file__), "..", "..", "scripts", "par")

SQ = [
    ("chr21", 46709983), ("chr22", 50818468), ("chrX", 156040895),
    ("chrY", 57227415), ("chrM", 16569),
]
UR = "UR:file:/ref/GRCh38.fasta"


def sq_lines():
    return ["@SQ\tSN:%s\tLN:%d\tM5:md5_%s\t%s" % (n, ln, n, UR) for n, ln in SQ]


# (contig, start, end); 1-based closed. Deliberately not grouped by contig order.
INTERVALS = [
    ("chr21", 10001, 11000),
    ("chr22", 20001, 21000),
    ("chrX", 10001, 11000),            # PAR1
    ("chrX", 2780001, 2781000),        # PAR1 (last)
    ("chrX", 2781001, 2782000),        # straddles PAR1 end -> dropped
    ("chrX", 2782001, 2783000),        # X-specific (first)
    ("chrX", 5000001, 5001000),        # X-specific
    ("chrX", 155700001, 155701000),    # X-specific (last)
    ("chrX", 155701001, 155702000),    # straddles PAR2 start -> dropped
    ("chrX", 155702001, 155703000),    # PAR2 (first)
    ("chrX", 156029001, 156030000),    # PAR2
    ("chrX", 156030001, 156031000),    # straddles PAR2 end -> dropped
    ("chrY", 10001, 11000),            # chrY PAR1 -> dropped
    ("chrY", 2781001, 2782000),        # straddles chrY PAR1 -> dropped
    ("chrY", 3000001, 3001000),        # kept
    ("chrY", 57000001, 57001000),      # chrY PAR2 -> dropped
    ("chrY", 57218001, 57219000),      # kept
    ("chrM", 1, 1000),
]
KEPT_X_SPECIFIC = [(2782001, 2783000), (5000001, 5001000), (155700001, 155701000)]


def write(path, text):
    with open(str(path), "w", newline="\n") as fh:
        fh.write(text)
    return str(path)


def read_lines(path):
    with open(str(path), newline="") as fh:
        return fh.read().splitlines()


def counts_file(path, intervals=INTERVALS):
    lines = ["@HD\tVN:1.6", "@RG\tID:GATKCopyNumber\tSM:sampleA"] + sq_lines()
    lines.append("CONTIG\tSTART\tEND\tCOUNT")
    lines += ["%s\t%d\t%d\t%d" % (c, s, e, i) for i, (c, s, e) in enumerate(intervals)]
    return write(path, "\n".join(lines) + "\n")


def interval_list_file(path, intervals=INTERVALS):
    lines = ["@HD\tVN:1.6"] + sq_lines()
    lines += ["%s\t%d\t%d\t+\t." % iv for iv in intervals]
    return write(path, "\n".join(lines) + "\n")


def annotated_file(path, intervals=INTERVALS):
    lines = ["@HD\tVN:1.6"] + sq_lines()
    lines.append("CONTIG\tSTART\tEND\tGC_CONTENT\tMAPPABILITY")
    lines += ["%s\t%d\t%d\t0.%06d\t1.000000" % (c, s, e, i) for i, (c, s, e) in enumerate(intervals)]
    return write(path, "\n".join(lines) + "\n")


def split_rows(lines):
    return [l.split("\t") for l in lines if not l.startswith("@") and not l.startswith("CONTIG")]


def relabel_to(tmp_path, builder, name):
    src = builder(tmp_path / (name + ".in"))
    dst = str(tmp_path / (name + ".out"))
    assert par_relabel.main(["--input", src, "--output", dst]) == 0
    return read_lines(dst)


# ------------------------------------------------------------------ dictionary

def test_pseudo_contigs_inserted_between_chr22_and_chrx():
    new = pc.pseudo_sq_lines(sq_lines())
    names = [pc.parse_sq(l)["SN"] for l in new]
    assert names == ["chr21", "chr22", "chrPAR1", "chrPAR2", "chrX", "chrY", "chrM"]
    par1 = pc.parse_sq(new[2])
    assert par1["LN"] == str(156040895)            # chrX length, no coordinate shift
    assert "M5" not in par1 and par1["UR"] == UR.split(":", 1)[1]
    assert pc.parse_sq(new[0])["M5"] == "md5_chr21"  # originals untouched
    assert pc.original_sq_lines(new) == sq_lines()


def test_pseudo_dict_refuses_already_relabeled():
    with pytest.raises(ValueError):
        pc.pseudo_sq_lines(pc.pseudo_sq_lines(sq_lines()))


def test_make_dict_cli(tmp_path):
    src = interval_list_file(tmp_path / "in.interval_list")
    out = str(tmp_path / "pseudo.dict")
    par_make_dict.main(["--input", src, "--output", out])
    lines = read_lines(out)
    assert lines[0] == "@HD\tVN:1.6"
    assert [pc.parse_sq(l)["SN"] for l in lines[1:]] == [
        "chr21", "chr22", "chrPAR1", "chrPAR2", "chrX", "chrY", "chrM"]


# ------------------------------------------------------------------ relabel

def test_counts_relabel(tmp_path):
    out = relabel_to(tmp_path, counts_file, "counts")
    # header preserved apart from @SQ
    assert out[0] == "@HD\tVN:1.6"
    assert out[1] == "@RG\tID:GATKCopyNumber\tSM:sampleA"
    sq = [l for l in out if l.startswith("@SQ")]
    assert [pc.parse_sq(l)["SN"] for l in sq] == [
        "chr21", "chr22", "chrPAR1", "chrPAR2", "chrX", "chrY", "chrM"]
    assert "CONTIG\tSTART\tEND\tCOUNT" in out
    rows = split_rows(out)
    contigs = [r[0] for r in rows]
    assert contigs == (["chr21", "chr22", "chrPAR1", "chrPAR1", "chrPAR2", "chrPAR2"]
                       + ["chrX"] * 3 + ["chrY", "chrY", "chrM"])
    by_contig = {}
    for r in rows:
        by_contig.setdefault(r[0], []).append((int(r[1]), int(r[2])))
    assert by_contig["chrPAR1"] == [(10001, 11000), (2780001, 2781000)]      # coordinates unchanged
    assert by_contig["chrPAR2"] == [(155702001, 155703000), (156029001, 156030000)]
    assert by_contig["chrX"] == KEPT_X_SPECIFIC                              # straddlers dropped
    assert by_contig["chrY"] == [(3000001, 3001000), (57218001, 57219000)]   # chrY PAR dropped
    # COUNT column untouched (row index of the original interval)
    first_par1 = [r for r in rows if r[0] == "chrPAR1"][0]
    assert first_par1[3] == str(INTERVALS.index(("chrX", 10001, 11000)))


def test_interval_list_relabel_has_no_column_header(tmp_path):
    out = relabel_to(tmp_path, interval_list_file, "il")
    assert not any(l.startswith("CONTIG") for l in out)
    rows = split_rows(out)
    assert rows[2] == ["chrPAR1", "10001", "11000", "+", "."]
    assert [r[0] for r in rows].count("chrPAR2") == 2


def test_annotated_relabel_keeps_annotations(tmp_path):
    out = relabel_to(tmp_path, annotated_file, "ann")
    assert "CONTIG\tSTART\tEND\tGC_CONTENT\tMAPPABILITY" in out
    rows = split_rows(out)
    par1 = [r for r in rows if r[0] == "chrPAR1"]
    i = INTERVALS.index(("chrX", 10001, 11000))
    assert par1[0] == ["chrPAR1", "10001", "11000", "0.%06d" % i, "1.000000"]


def test_three_formats_make_identical_decisions(tmp_path):
    results = []
    for builder, name in ((counts_file, "c"), (interval_list_file, "i"), (annotated_file, "a")):
        rows = split_rows(relabel_to(tmp_path, builder, name))
        results.append([tuple(r[:3]) for r in rows])
    assert results[0] == results[1] == results[2]


def test_output_is_sorted_in_new_dictionary_order_even_for_unsorted_input(tmp_path):
    shuffled = list(reversed(INTERVALS))
    src = counts_file(tmp_path / "in.tsv", shuffled)
    dst = str(tmp_path / "out.tsv")
    par_relabel.main(["--input", src, "--output", dst])
    rows = split_rows(read_lines(dst))
    order = {n: i for i, n in enumerate(["chr21", "chr22", "chrPAR1", "chrPAR2", "chrX", "chrY", "chrM"])}
    keys = [(order[r[0]], int(r[1])) for r in rows]
    assert keys == sorted(keys)


def test_exclude_mode_drops_all_par_and_keeps_dictionary(tmp_path):
    src = counts_file(tmp_path / "in.tsv")
    dst = str(tmp_path / "out.tsv")
    par_relabel.main(["--input", src, "--output", dst, "--mode", "exclude"])
    out = read_lines(dst)
    assert [l for l in out if l.startswith("@SQ")] == sq_lines()
    rows = split_rows(out)
    assert [(r[0], int(r[1])) for r in rows if r[0] == "chrX"] == [("chrX", s) for s, _ in KEPT_X_SPECIFIC]
    assert [r[0] for r in rows].count("chrY") == 2
    assert not any(r[0].startswith("chrPAR") for r in rows)


def test_exclude_bed(tmp_path):
    bed = write(tmp_path / "xtr.bed", "# comment\nchrX\t4999999\t5000500\nchrX\t5000400\t5000600\n")
    src = interval_list_file(tmp_path / "in.il")
    dst = str(tmp_path / "out.il")
    par_relabel.main(["--input", src, "--output", dst, "--exclude-bed", bed])
    rows = split_rows(read_lines(dst))
    x = [(int(r[1]), int(r[2])) for r in rows if r[0] == "chrX"]
    assert (5000001, 5001000) not in x and (2782001, 2783000) in x


def test_bed_overlap_edges():
    bed = {"chrX": ([100], [200])}                     # 0-based [100,200) = 1-based 101..200
    assert pc.overlaps_bed(bed, "chrX", 200, 300)
    assert not pc.overlaps_bed(bed, "chrX", 201, 300)
    assert pc.overlaps_bed(bed, "chrX", 1, 101)
    assert not pc.overlaps_bed(bed, "chrX", 1, 100)
    assert not pc.overlaps_bed(bed, "chr1", 1, 1000)


def test_classify_boundaries():
    c = pc.classify_interval
    assert c("chrX", 10001, 2781479, "pseudo_contig") == ("chrPAR1", "par1")
    assert c("chrX", 10001, 2781480, "pseudo_contig")[0] is None
    assert c("chrX", 155701383, 156030895, "pseudo_contig") == ("chrPAR2", "par2")
    assert c("chrY", 56887903, 56887903, "pseudo_contig")[0] is None
    assert c("chr1", 10001, 20000, "pseudo_contig") == ("chr1", "kept")


# ------------------------------------------------------------------ reconstruct

PSEUDO_ORDER = ["chr21", "chr22", "chrPAR1", "chrPAR2", "chrX", "chrY", "chrM"]


def vcf_text(records, pseudo=True):
    meta = ["##fileformat=VCFv4.2",
            "##ALT=<ID=DEL,Description=\"Deletion\">",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">"]
    names = PSEUDO_ORDER if pseudo else [n for n, _ in SQ]
    lengths = dict(SQ)
    lengths["chrPAR1"] = lengths["chrPAR2"] = lengths["chrX"]
    meta += ["##contig=<ID=%s,length=%d>" % (n, lengths[n]) for n in names]
    meta.append("##GATKCommandLine=<ID=Postprocess>")
    cols = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsampleA"
    lines = meta + [cols]
    for contig, pos, end, alt, gt in records:
        lines.append("%s\t%d\t.\tN\t%s\t.\t.\tEND=%d\tGT:CN\t%s:2" % (contig, pos, alt, end, gt))
    return "\n".join(lines) + "\n"


def write_gz(path, text):
    with gzip.open(str(path), "wt", newline="\n") as fh:
        fh.write(text)
    return str(path)


# fake postprocess output, written (like GATK) in pseudo dictionary order
FAKE_INTERVALS = [
    ("chr21", 10001, 11000, "<DEL>", "0"),
    ("chrPAR1", 10001, 11000, ".", "0"),
    ("chrPAR1", 2780001, 2781000, "<DEL>", "1"),
    ("chrPAR2", 155702001, 155703000, "<DEL>", "1"),
    ("chrPAR2", 156029001, 156030000, ".", "0"),
    ("chrX", 2782001, 2783000, "<DEL>", "1"),
    ("chrX", 5000001, 5001000, ".", "0"),
    ("chrX", 155700001, 155701000, "<DEL>", "1"),
]
FAKE_SEGMENTS = [
    ("chr21", 10001, 11000, ".", "0"),
    ("chrPAR1", 10001, 2779000, ".", "0"),                 # reference
    ("chrPAR1", 2780001, 2781000, "<DEL>", "1"),           # non-ref, ends at last PAR1 interval
    ("chrPAR2", 155702001, 155703000, "<DEL>", "1"),       # non-ref, starts at first PAR2 interval
    ("chrPAR2", 156029001, 156030000, ".", "0"),
    ("chrX", 2782001, 2783000, "<DEL>", "1"),              # non-ref, first X-specific interval
    ("chrX", 5000001, 5001000, ".", "0"),
    ("chrX", 155700001, 155701000, "<DEL>", "1"),          # non-ref, last X-specific interval
]


def run_reconstruct(tmp_path, segments=FAKE_SEGMENTS, intervals=FAKE_INTERVALS, ploidy_x=1):
    pseudo_dict = write(tmp_path / "pseudo.dict", "\n".join(
        ["@HD\tVN:1.6"] + pc.pseudo_sq_lines(sq_lines())) + "\n")
    ploidy = write(tmp_path / "contig_ploidy.tsv",
                   "@RG\tID:GATKCopyNumber\tSM:sampleA\nCONTIG\tPLOIDY\tPLOIDY_GQ\n"
                   "chr21\t2\t99\nchr22\t2\t99\nchrPAR1\t2\t99\nchrPAR2\t2\t99\n"
                   "chrX\t%d\t99\nchrY\t%d\t99\nchrM\t2\t99\n" % (ploidy_x, 2 - ploidy_x))
    seg_in = write_gz(tmp_path / "seg.vcf.gz", vcf_text(segments))
    int_in = write_gz(tmp_path / "int.vcf.gz", vcf_text(intervals))
    dcr_in = write(tmp_path / "dcr.tsv", "@RG\tID:GATKCopyNumber\tSM:sampleA\n"
                   + "\n".join(pc.pseudo_sq_lines(sq_lines())) + "\n"
                   "CONTIG\tSTART\tEND\tLOG2_COPY_RATIO\n"
                   "chr21\t10001\t11000\t0.0\nchrPAR1\t10001\t11000\t0.01\nchrPAR2\t155702001\t155703000\t-0.02\n"
                   "chrX\t5000001\t5001000\t-1.0\n")
    out = dict(seg=str(tmp_path / "seg.out.vcf"), ints=str(tmp_path / "int.out.vcf"),
               dcr=str(tmp_path / "dcr.out.tsv"))
    par_reconstruct.main(["--pseudo-dict", pseudo_dict, "--contig-ploidy", ploidy,
                          "--segments-vcf", seg_in, "--out-segments-vcf", out["seg"],
                          "--intervals-vcf", int_in, "--out-intervals-vcf", out["ints"],
                          "--denoised-copy-ratios", dcr_in, "--out-denoised-copy-ratios", out["dcr"]])
    return out


def vcf_records(path):
    return [l.split("\t") for l in read_lines(path) if not l.startswith("#")]


def info_of(rec):
    return rec[7]


def test_round_trip_coordinates_and_header(tmp_path):
    # relabel the original intervals, then pretend GATK produced one record per kept interval
    src = interval_list_file(tmp_path / "in.il")
    dst = str(tmp_path / "relabeled.il")
    par_relabel.main(["--input", src, "--output", dst])
    kept = [(r[0], int(r[1]), int(r[2]), "<DEL>", "1") for r in split_rows(read_lines(dst))]
    out = run_reconstruct(tmp_path, segments=kept, intervals=kept)

    originals = {(c, s, e) for c, s, e in INTERVALS if c == "chrX"}
    rebuilt = {(r[0], int(r[1]), int(r[7].split("END=")[1].split(";")[0])) for r in vcf_records(out["ints"])
               if r[0] == "chrX"}
    # every reconstructed chrX record equals an original chrX interval; dropped straddlers absent
    assert rebuilt <= originals
    assert len(rebuilt) == 3 + 2 + 2          # X-specific + PAR1 + PAR2 intervals
    assert ("chrX", 2781001, 2782000) not in rebuilt

    lines = read_lines(out["ints"])
    assert not any("chrPAR" in l for l in lines)
    contig_ids = [l.split("ID=")[1].split(",")[0] for l in lines if l.startswith("##contig=")]
    assert contig_ids == [n for n, _ in SQ]
    # contig block sits where it was; other header lines preserved
    assert lines[0] == "##fileformat=VCFv4.2"
    assert "##GATKCommandLine=<ID=Postprocess>" in lines


def test_vcf_sorted_in_original_order(tmp_path):
    out = run_reconstruct(tmp_path)
    recs = vcf_records(out["seg"])
    assert [r[0] for r in recs] == ["chr21"] + ["chrX"] * 7
    xpos = [int(r[1]) for r in recs if r[0] == "chrX"]
    assert xpos == sorted(xpos)
    assert xpos[0] == 10001 and xpos[-1] == 156029001     # PAR1 first, PAR2 last on chrX


def test_info_tags_par_and_refcn(tmp_path):
    out = run_reconstruct(tmp_path, ploidy_x=1)
    recs = {(r[0], int(r[1])): info_of(r) for r in vcf_records(out["ints"])}
    assert recs[("chrX", 10001)] == "END=11000;PAR=PAR1;REFCN=2"
    assert recs[("chrX", 155702001)].endswith("PAR=PAR2;REFCN=2")
    assert recs[("chrX", 5000001)] == "END=5001000;REFCN=1"     # chrX ploidy of the sample
    assert recs[("chr21", 10001)] == "END=11000"                # untouched
    header = read_lines(out["ints"])
    assert any(l.startswith("##INFO=<ID=PAR,") for l in header)
    assert any(l.startswith("##INFO=<ID=REFCN,") for l in header)
    assert not any("PAR_BOUNDARY_ADJACENT" in l for l in header)  # segments only


def test_refcn_follows_female_ploidy(tmp_path):
    out = run_reconstruct(tmp_path, ploidy_x=2)
    recs = {(r[0], int(r[1])): info_of(r) for r in vcf_records(out["seg"])}
    assert "REFCN=2" in recs[("chrX", 5000001)] and "PAR" not in recs[("chrX", 5000001)]


def test_boundary_adjacent_flag(tmp_path):
    out = run_reconstruct(tmp_path)
    recs = {(r[0], int(r[1])): info_of(r) for r in vcf_records(out["seg"])}
    flagged = {k for k, v in recs.items() if "PAR_BOUNDARY_ADJACENT" in v}
    assert flagged == {
        ("chrX", 2780001),     # PAR1 side: non-ref, ends at last PAR1 interval
        ("chrX", 2782001),     # X side after PAR1: non-ref, starts at first X-specific interval
        ("chrX", 155700001),   # X side before PAR2: non-ref, ends at last X-specific interval
        ("chrX", 155702001),   # PAR2 side: non-ref, starts at first PAR2 interval
    }
    # reference segments and mid-PAR / mid-X records are not flagged
    assert "PAR_BOUNDARY_ADJACENT" not in recs[("chrX", 10001)]
    assert "PAR_BOUNDARY_ADJACENT" not in recs[("chrX", 5000001)]
    assert any(l.startswith("##INFO=<ID=PAR_BOUNDARY_ADJACENT,") for l in read_lines(out["seg"]))


def test_boundary_flag_requires_non_reference(tmp_path):
    segs = [("chrPAR1", 2780001, 2781000, ".", "0"), ("chrX", 2782001, 2783000, ".", "0")]
    out = run_reconstruct(tmp_path, segments=segs)
    assert not any("PAR_BOUNDARY_ADJACENT" in info_of(r) for r in vcf_records(out["seg"]))


def test_boundary_uses_present_intervals_not_nominal_boundary(tmp_path):
    # PAR1 interval at 2780001-2781000 exists, a segment that stops one interval earlier is no flag
    segs = [("chrPAR1", 10001, 11000, "<DEL>", "1")]
    out = run_reconstruct(tmp_path, segments=segs)
    assert "PAR_BOUNDARY_ADJACENT" not in info_of(vcf_records(out["seg"])[0])


def test_denoised_copy_ratios_reconstruct(tmp_path):
    out = run_reconstruct(tmp_path)
    lines = read_lines(out["dcr"])
    assert [l for l in lines if l.startswith("@SQ")] == sq_lines()
    assert lines[0] == "@RG\tID:GATKCopyNumber\tSM:sampleA"
    rows = split_rows(lines)
    assert [(r[0], r[1]) for r in rows] == [
        ("chr21", "10001"), ("chrX", "10001"), ("chrX", "5000001"), ("chrX", "155702001")]
    assert rows[1][3] == "0.01"


def test_non_reference_logic():
    f = lambda alt, gt: ["c", "1", ".", "N", alt, ".", ".", "END=2", "GT:CN", gt + ":2"]
    assert not pc.is_non_reference(f(".", "0"))
    assert not pc.is_non_reference(f("<DEL>", "0/0"))
    assert pc.is_non_reference(f("<DEL>", "0/1"))
    assert pc.is_non_reference(f("<DUP>", "1"))
    assert pc.is_non_reference(f("<DEL>", "./."))     # as kept by the finalize rule
    assert not pc.is_non_reference(f(".", "1"))


# ------------------------------------------------------------------ ploidy QC

def ploidy_file(path, x, y, par1, par2=None):
    par2 = par1 if par2 is None else par2
    return write(path, "@RG\tID:GATKCopyNumber\tSM:S1\nCONTIG\tPLOIDY\tPLOIDY_GQ\n"
                 "chr21\t2\t99.0\nchrPAR1\t%d\t80.5\nchrPAR2\t%d\t70.0\nchrX\t%d\t99.0\nchrY\t%d\t60.0\n"
                 % (par1, par2, x, y))


def qc(tmp_path, x, y, par1, par2=None):
    src = ploidy_file(tmp_path / "contig_ploidy.tsv", x, y, par1, par2)
    dst = str(tmp_path / "qc.tsv")
    assert par_ploidy_qc.main(["--contig-ploidy", src, "--output", dst]) == 0
    lines = read_lines(dst)
    assert lines[0].split("\t") == list(pc.QC_HEADER)
    return [l.split("\t") for l in lines[1:]]


def test_qc_46xx_and_46xy_ok(tmp_path):
    rows = qc(tmp_path, 2, 0, 2)
    assert rows[0] == ["S1", "chrPAR1", "2", "2", "80.5", "OK"]
    assert [r[-1] for r in qc(tmp_path, 1, 1, 2)] == ["OK", "OK"]


def test_qc_45x_ok(tmp_path):
    rows = qc(tmp_path, 1, 0, 1)
    assert [r[-1] for r in rows] == ["OK", "OK"]


def test_qc_mismatch_does_not_fail(tmp_path, capsys):
    rows = qc(tmp_path, 1, 1, 1, 2)            # 46,XY with PAR1 called as ploidy 1
    assert rows[0][1:] == ["chrPAR1", "1", "2", "80.5", "MISMATCH"]
    assert rows[1][-1] == "OK"
    assert "MISMATCH" in capsys.readouterr().err


def test_qc_na_when_contigs_missing(tmp_path):
    src = write(tmp_path / "p.tsv", "@RG\tID:GATKCopyNumber\tSM:S1\nCONTIG\tPLOIDY\tPLOIDY_GQ\nchr21\t2\t99\n")
    dst = str(tmp_path / "qc.tsv")
    par_ploidy_qc.main(["--contig-ploidy", src, "--output", dst])
    assert all(l.split("\t")[-1] == "NA" for l in read_lines(dst)[1:])


# ------------------------------------------------------------------ CLI smoke

@pytest.mark.parametrize("script", ["par_relabel.py", "par_make_dict.py", "par_reconstruct.py",
                                    "par_ploidy_qc.py"])
def test_help(script):
    res = subprocess.run([sys.executable, os.path.join(SCRIPTS, script), "--help"],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    assert res.returncode == 0 and "usage" in res.stdout.lower()
