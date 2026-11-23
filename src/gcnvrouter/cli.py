"""Command line: `gcnvr explore ...` (or `python -m gcnvrouter.cli explore ...`)."""

import argparse
import json
import sys

from .explore import explore, write_outputs
from .io_counts import find_count_files
from .profiles import PROFILES, get_profile


def _explore(args) -> int:
    if args.counts_list:
        with open(args.counts_list) as fh:
            paths = [l.strip() for l in fh if l.strip() and not l.startswith("#")]
    else:
        paths = find_count_files(args.counts_dir, args.pattern)
    profile = get_profile(args.profile).with_overrides(explore_n_min=args.n_min, k_pc_max=args.k_pc_max)
    exp = explore(paths, profile, metadata_path=args.metadata, annotated_path=args.annotated,
                  mask_bed_paths=args.mask_bed or (), exclude=args.exclude or (), seed=args.seed,
                  n_boot=args.n_boot)
    summ = write_outputs(exp, args.output, html=not args.no_html)
    print(json.dumps({k: summ[k] for k in ("n_samples", "n_fit", "preflight", "n_features",
                                           "k_parallel_analysis", "K", "covariance", "cluster_sizes",
                                           "bootstrap_median_ari", "status_counts", "warnings")}, indent=2))
    print("outputs written to %s" % args.output)
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="gcnvr", description="Batch discovery and routing for GATK gCNV")
    sub = ap.add_subparsers(dest="cmd", required=True)

    e = sub.add_parser("explore", help="cluster historical samples and write a report")
    src = e.add_mutually_exclusive_group(required=True)
    src.add_argument("--counts-dir", help="directory of CollectReadCounts TSV files")
    src.add_argument("--counts-list", help="file with one count-file path per line")
    e.add_argument("--pattern", default="*.tsv", help="glob within --counts-dir (default *.tsv)")
    e.add_argument("--profile", required=True, choices=sorted(PROFILES))
    e.add_argument("-o", "--output", required=True, help="output directory")
    e.add_argument("--metadata", help="TSV with sample_id + metadata columns (run, kit, date, ...)")
    e.add_argument("--annotated", help="GATK AnnotateIntervals TSV (GC/mappability/segdup)")
    e.add_argument("--mask-bed", action="append", help="BED of regions to exclude from features (repeatable)")
    e.add_argument("--exclude", action="append", help="sample to leave out of the fit (repeatable)")
    e.add_argument("--n-min", type=int, help="minimum cluster size during exploration (profile default)")
    e.add_argument("--k-pc-max", type=int, help="maximum number of PCs (profile default)")
    e.add_argument("--seed", type=int, default=1729)
    e.add_argument("--n-boot", type=int, default=100)
    e.add_argument("--no-html", action="store_true")
    e.set_defaults(func=_explore)

    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
