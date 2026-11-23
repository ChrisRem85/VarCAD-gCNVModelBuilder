#!/usr/bin/env python3
"""Relabel PAR1/PAR2 of chrX as pseudo-contigs chrPAR1/chrPAR2 (or exclude PAR intervals).

One implementation for all three GATK interval-style formats (detected from the file):
  * CollectReadCounts TSV      '@' header, 'CONTIG START END COUNT', rows
  * annotated-intervals TSV    '@' header, 'CONTIG START END <annotations>', rows
  * Picard interval_list       '@' header, rows 'contig start end strand name'

Modes
  pseudo_contig  chrX intervals fully inside PAR1/PAR2 -> chrPAR1/chrPAR2 (coordinates
                 unchanged); chrX intervals straddling a PAR boundary are dropped; all chrY
                 intervals overlapping a chrY PAR are dropped. The @SQ header is replaced by
                 the dictionary with chrPAR1/chrPAR2 inserted before chrX, and the rows are
                 re-sorted in that dictionary order.
  exclude        every chrX/chrY interval overlapping a PAR is dropped; header and contig
                 order are unchanged.

--exclude-bed additionally drops intervals overlapping a BED region (e.g. XTR, 0-based BED
coordinates, original contig names). No coordinates are hard-coded for it.
"""

import argparse
import sys
from collections import Counter

import par_common as pc


def relabel(header, column_header, rows, mode, exclude_bed=None):
    """Return (new_header, new_rows, stats)."""
    stats = Counter()
    sq_lines = pc.sq_lines_of(header)
    if mode == pc.MODE_PSEUDO:
        if not sq_lines:
            raise ValueError("input has no @SQ header lines; cannot build the pseudo dictionary")
        new_sq = pc.pseudo_sq_lines(sq_lines)
        new_header = pc.replace_sq_block(header, new_sq)
        order = pc.contig_order(new_sq)
    else:
        new_header = list(header)
        order = None

    buckets = {}   # pseudo mode: contig index -> [(start, line)]
    kept_rows = []  # exclude mode: input order is kept
    for line in rows:
        fields = line.split("\t", 3)
        contig, start, end = fields[0], int(fields[1]), int(fields[2])
        if exclude_bed and pc.overlaps_bed(exclude_bed, contig, start, end):
            stats["dropped_exclude_bed"] += 1
            continue
        new_contig, reason = pc.classify_interval(contig, start, end, mode)
        stats[reason] += 1
        if new_contig is None:
            continue
        if new_contig != contig:
            line = new_contig + line[len(contig):]
        if order is None:
            kept_rows.append(line)
        else:
            if new_contig not in order:
                raise ValueError("contig %r is not in the @SQ dictionary" % new_contig)
            buckets.setdefault(order[new_contig], []).append((start, line))

    if order is None:
        return new_header, kept_rows, stats
    new_rows = []
    for idx in sorted(buckets):
        entries = buckets[idx]
        entries.sort(key=lambda e: e[0])  # stable; already sorted for GATK-produced input
        new_rows.extend(line for _, line in entries)
    return new_header, new_rows, stats


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="counts TSV, annotated TSV or interval_list")
    ap.add_argument("--output", required=True)
    ap.add_argument("--mode", choices=[pc.MODE_PSEUDO, pc.MODE_EXCLUDE], default=pc.MODE_PSEUDO)
    ap.add_argument("--exclude-bed", default=None,
                    help="optional BED of regions whose intervals are dropped (e.g. XTR)")
    args = ap.parse_args(argv)

    header, column_header, rows = pc.load_table(args.input)
    bed = pc.load_bed(args.exclude_bed) if args.exclude_bed else None
    new_header, new_rows, stats = relabel(header, column_header, rows, args.mode, bed)
    pc.write_table(args.output, new_header, column_header, new_rows)

    sys.stderr.write("par_relabel %s: %d rows in, %d rows out; %s\n" % (
        args.input, len(rows), len(new_rows),
        ", ".join("%s=%d" % kv for kv in sorted(stats.items()))))
    return 0


if __name__ == "__main__":
    sys.exit(main())
