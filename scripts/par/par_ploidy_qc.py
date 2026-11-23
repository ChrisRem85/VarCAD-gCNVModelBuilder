#!/usr/bin/env python3
"""PAR ploidy consistency QC for one sample.

Reads the sample's contig_ploidy.tsv (DetermineGermlineContigPloidy output) and checks

    ploidy(chrPAR1) == ploidy(chrX) + ploidy(chrY)
    ploidy(chrPAR2) == ploidy(chrX) + ploidy(chrY)

Writes a TSV with columns: sample, contig, observed, expected, PLOIDY_GQ, status.
status is OK, MISMATCH (whole-PAR copy-number change, sample swap, sex-check problem), or NA
when chrX/chrY/chrPARn are missing from the ploidy calls. A mismatch is a finding, not an
error: the exit code is 0 either way.
"""

import argparse
import sys

import par_common as pc


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--contig-ploidy", required=True, help="contig_ploidy.tsv of the sample")
    ap.add_argument("--output", required=True, help="QC TSV to write")
    args = ap.parse_args(argv)

    sample, ploidy = pc.read_contig_ploidy(args.contig_ploidy)
    rows = pc.par_ploidy_qc_rows(sample, ploidy)
    with pc.open_text(args.output, "w") as out:
        out.write("\t".join(pc.QC_HEADER) + "\n")
        for row in rows:
            out.write("\t".join(row) + "\n")
    for row in rows:
        if row[-1] == "MISMATCH":
            sys.stderr.write("par_ploidy_qc: MISMATCH sample=%s %s observed=%s expected=%s\n"
                             % (row[0], row[1], row[2], row[3]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
