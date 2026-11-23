"""Shared helpers for the PAR1/PAR2 pseudo-contig workflow.

Python 3.6+, standard library only (the scripts run inside the broadinstitute/gatk
container).

Coordinates are GRCh38, 1-based, closed.
"""

import bisect
import gzip
import re
import sys
from collections import OrderedDict

# --------------------------------------------------------------------------------------
# Coordinates
# --------------------------------------------------------------------------------------

PAR_COORDS = {
    "chrX": OrderedDict([("PAR1", (10001, 2781479)), ("PAR2", (155701383, 156030895))]),
    "chrY": OrderedDict([("PAR1", (10001, 2781479)), ("PAR2", (56887903, 57217415))]),
}

PSEUDO_CONTIGS = ("chrPAR1", "chrPAR2")
PSEUDO_FOR_REGION = {"PAR1": "chrPAR1", "PAR2": "chrPAR2"}
REGION_FOR_PSEUDO = {"chrPAR1": "PAR1", "chrPAR2": "PAR2"}

MODE_PSEUDO = "pseudo_contig"
MODE_EXCLUDE = "exclude"


# --------------------------------------------------------------------------------------
# File helpers
# --------------------------------------------------------------------------------------

def open_text(path, mode="r"):
    """Open a text file for reading/writing. '-' is stdin/stdout, '*.gz' is gzip."""
    if path == "-":
        return sys.stdin if mode == "r" else sys.stdout
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t", encoding="utf-8", newline="" if mode == "r" else "\n")
    if mode == "r":
        return open(path, "r", encoding="utf-8")
    return open(path, "w", encoding="utf-8", newline="\n")


def load_table(path):
    """Read a GATK-style table: '@' header lines, optional 'CONTIG ...' column header, rows.

    Handles CollectReadCounts TSV, annotated-intervals TSV (column header present) and
    Picard interval_list (no column header).
    Returns (header_lines, column_header_or_None, row_lines); line terminators removed.
    """
    header, column_header, rows = [], None, []
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith("@"):
                header.append(line)
            elif column_header is None and not rows and line.split("\t", 1)[0] == "CONTIG":
                column_header = line
            else:
                rows.append(line)
    return header, column_header, rows


def write_table(path, header, column_header, rows):
    with open_text(path, "w") as out:
        for line in header:
            out.write(line + "\n")
        if column_header is not None:
            out.write(column_header + "\n")
        for line in rows:
            out.write(line + "\n")


# --------------------------------------------------------------------------------------
# Sequence dictionary
# --------------------------------------------------------------------------------------

def parse_sq(line):
    """'@SQ\\tSN:chr1\\tLN:10' -> OrderedDict(SN=chr1, LN=10)."""
    fields = line.rstrip("\r\n").split("\t")
    if fields[0] != "@SQ":
        raise ValueError("not an @SQ line: %r" % line)
    tags = OrderedDict()
    for field in fields[1:]:
        key, _, value = field.partition(":")
        tags[key] = value
    return tags


def format_sq(tags):
    return "@SQ\t" + "\t".join("%s:%s" % (k, v) for k, v in tags.items())


def sq_lines_of(header):
    return [line for line in header if line.startswith("@SQ\t")]


def pseudo_sq_lines(sq_lines):
    """Insert chrPAR1/chrPAR2 immediately before chrX (= right after chr22 in GRCh38).

    Both get LN of chrX; other tags (e.g. UR) are copied from chrX, M5 is dropped because
    the checksum of chrX does not describe the pseudo-contig.
    """
    names = [parse_sq(line)["SN"] for line in sq_lines]
    if any(name in PSEUDO_CONTIGS for name in names):
        raise ValueError("dictionary already contains chrPAR1/chrPAR2 (input already relabeled?)")
    if "chrX" not in names:
        raise ValueError("dictionary has no chrX contig")
    x_index = names.index("chrX")
    x_tags = parse_sq(sq_lines[x_index])
    new_lines = []
    for contig in PSEUDO_CONTIGS:
        tags = OrderedDict(x_tags)
        tags["SN"] = contig
        tags.pop("M5", None)
        new_lines.append(format_sq(tags))
    return list(sq_lines[:x_index]) + new_lines + list(sq_lines[x_index:])


def original_sq_lines(sq_lines):
    """Remove chrPAR1/chrPAR2 again."""
    return [line for line in sq_lines if parse_sq(line)["SN"] not in PSEUDO_CONTIGS]


def contig_order(sq_lines):
    return dict((parse_sq(line)["SN"], i) for i, line in enumerate(sq_lines))


def contig_lengths(sq_lines):
    out = OrderedDict()
    for line in sq_lines:
        tags = parse_sq(line)
        out[tags["SN"]] = int(tags["LN"])
    return out


def replace_sq_block(header, new_sq_lines):
    """Replace all @SQ lines by new_sq_lines at the position of the first @SQ line."""
    out, done = [], False
    for line in header:
        if line.startswith("@SQ\t"):
            if not done:
                out.extend(new_sq_lines)
                done = True
        else:
            out.append(line)
    if not done:
        out.extend(new_sq_lines)
    return out


def read_header(path):
    """Leading '@' lines of a .dict / interval_list / counts / annotated file."""
    header = []
    with open_text(path) as fh:
        for line in fh:
            if not line.startswith("@"):
                break
            header.append(line.rstrip("\r\n"))
    return header


# --------------------------------------------------------------------------------------
# Interval classification
# --------------------------------------------------------------------------------------

def _overlaps(start, end, region):
    return start <= region[1] and end >= region[0]


def _inside(start, end, region):
    return start >= region[0] and end <= region[1]


def classify_interval(contig, start, end, mode):
    """Decide what happens to one 1-based closed interval.

    Returns (new_contig, reason); new_contig is None when the interval is dropped.
    Reasons: kept, par1, par2, straddle, chrY_par, par_excluded.
    """
    if contig not in PAR_COORDS:
        return contig, "kept"
    regions = PAR_COORDS[contig]
    if mode == MODE_EXCLUDE:
        for region in regions.values():
            if _overlaps(start, end, region):
                return None, "par_excluded"
        return contig, "kept"
    if mode != MODE_PSEUDO:
        raise ValueError("unknown PAR mode: %s" % mode)
    if contig == "chrY":
        for region in regions.values():
            if _overlaps(start, end, region):
                return None, "chrY_par"
        return contig, "kept"
    for name, region in regions.items():  # chrX
        if _inside(start, end, region):
            return PSEUDO_FOR_REGION[name], name.lower()
        if _overlaps(start, end, region):
            return None, "straddle"
    return contig, "kept"


# --------------------------------------------------------------------------------------
# Optional exclusion BED (0-based, half-open on disk)
# --------------------------------------------------------------------------------------

def load_bed(path):
    """Return {contig: (starts, ends)} with merged, sorted regions (0-based half-open)."""
    raw = {}
    with open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split()
            raw.setdefault(fields[0], []).append((int(fields[1]), int(fields[2])))
    merged = {}
    for contig, regions in raw.items():
        regions.sort()
        out = []
        for start, end in regions:
            if out and start <= out[-1][1]:
                out[-1] = (out[-1][0], max(out[-1][1], end))
            else:
                out.append((start, end))
        merged[contig] = ([r[0] for r in out], [r[1] for r in out])
    return merged


def overlaps_bed(bed, contig, start, end):
    """Does the 1-based closed interval [start, end] overlap any BED region?"""
    if contig not in bed:
        return False
    starts, ends = bed[contig]
    i = bisect.bisect_left(ends, start)  # first region with bed_end >= start
    return i < len(starts) and starts[i] < end


# --------------------------------------------------------------------------------------
# Contig ploidy calls (DetermineGermlineContigPloidy output)
# --------------------------------------------------------------------------------------

def read_contig_ploidy(path):
    """Return (sample_name, OrderedDict contig -> (ploidy:int, ploidy_gq:str))."""
    sample = None
    ploidy = OrderedDict()
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith("@RG"):
                for field in line.split("\t")[1:]:
                    if field.startswith("SM:"):
                        sample = field[3:]
                continue
            if line.startswith("@"):
                continue
            fields = line.split("\t")
            if fields[0] == "CONTIG":
                continue
            gq = fields[2] if len(fields) > 2 else "NA"
            ploidy[fields[0]] = (int(round(float(fields[1]))), gq)
    return (sample if sample is not None else "NA"), ploidy


def par_ploidy_qc_rows(sample, ploidy):
    """Expected ploidy(chrPARn) = ploidy(chrX) + ploidy(chrY).

    Returns rows (sample, contig, observed, expected, PLOIDY_GQ, status); status is
    OK, MISMATCH, or NA when a needed contig is missing from the ploidy calls.
    """
    rows = []
    for pseudo in PSEUDO_CONTIGS:
        observed = ploidy.get(pseudo)
        x, y = ploidy.get("chrX"), ploidy.get("chrY")
        if observed is None or x is None or y is None:
            rows.append((sample, pseudo,
                         "NA" if observed is None else str(observed[0]), "NA",
                         "NA" if observed is None else observed[1], "NA"))
            continue
        expected = x[0] + y[0]
        status = "OK" if observed[0] == expected else "MISMATCH"
        rows.append((sample, pseudo, str(observed[0]), str(expected), observed[1], status))
    return rows


QC_HEADER = ("sample", "contig", "observed", "expected", "PLOIDY_GQ", "status")


# --------------------------------------------------------------------------------------
# VCF helpers
# --------------------------------------------------------------------------------------

_END_RE = re.compile(r"(?:^|;)END=(\d+)")


def info_end(info, pos):
    match = _END_RE.search(info)
    return int(match.group(1)) if match else pos


def is_non_reference(fields):
    """Non-reference segment/interval record, consistent with finalize_genotyped_segments.

    That rule keeps records matching  GT != "ref" && ALT != "."  in bcftools terms: ALT must
    be present and the genotype must not be all-0 (a missing './.' genotype counts as
    non-reference there, so it does here).
    """
    if fields[4] == ".":
        return False
    if len(fields) < 10:
        return True
    keys = fields[8].split(":")
    if "GT" not in keys:
        return True
    values = fields[9].split(":")
    idx = keys.index("GT")
    if idx >= len(values):
        return True
    alleles = re.split(r"[/|]", values[idx])
    return not all(a == "0" for a in alleles)
