#!/usr/bin/env python3
"""Write the pseudo-contig sequence dictionary (SAM-header style .dict).

The dictionary is derived from the @SQ header of an existing file (preprocessed
interval_list, counts TSV, annotated intervals or a .dict), so it always matches the files
that par_relabel.py relabels: chrPAR1 and chrPAR2 are inserted immediately before chrX with
LN = chrX length. M5/UR (and other tags) of the original contigs are carried over; the new
contigs get no M5. Coordinates are not shifted.

Output: the input's @HD line (or '@HD VN:1.6') followed by the @SQ lines.
"""

import argparse
import sys

import par_common as pc


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="file whose @SQ header is the source")
    ap.add_argument("--output", required=True, help="pseudo .dict to write")
    args = ap.parse_args(argv)

    header = pc.read_header(args.input)
    sq_lines = pc.sq_lines_of(header)
    if not sq_lines:
        sys.exit("par_make_dict: no @SQ lines in %s" % args.input)
    hd = [line for line in header if line.startswith("@HD")] or ["@HD\tVN:1.6"]
    new_sq = pc.pseudo_sq_lines(sq_lines)
    with pc.open_text(args.output, "w") as out:
        for line in hd[:1] + new_sq:
            out.write(line + "\n")
    sys.stderr.write("par_make_dict: wrote %d @SQ lines to %s\n" % (len(new_sq), args.output))
    return 0


if __name__ == "__main__":
    sys.exit(main())
