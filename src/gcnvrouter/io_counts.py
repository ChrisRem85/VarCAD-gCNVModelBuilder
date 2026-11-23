"""Reading GATK CollectReadCounts TSV files into one aligned matrix."""

import glob
import os
from dataclasses import dataclass
from typing import List, Sequence, Tuple

import numpy as np
import pandas as pd

COORDS = ["CONTIG", "START", "END"]


@dataclass
class CountMatrix:
    samples: List[str]
    intervals: pd.DataFrame          # CONTIG, START, END (1-based, closed)
    counts: np.ndarray               # n_samples x n_intervals, float64

    def subset(self, mask) -> "CountMatrix":
        idx = np.flatnonzero(mask)
        return CountMatrix([self.samples[i] for i in idx], self.intervals, self.counts[idx])


def _sample_from_header(lines: Sequence[str]) -> str:
    for line in lines:
        if line.startswith("@RG"):
            for field in line.rstrip("\n").split("\t"):
                if field.startswith("SM:"):
                    return field[3:]
    return ""


def read_count_file(path: str) -> Tuple[str, pd.DataFrame, np.ndarray]:
    """Return (sample name, intervals, counts) for one CollectReadCounts TSV."""
    header = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("@"):
                header.append(line)
            else:
                break
    table = pd.read_csv(path, sep="\t", skiprows=len(header),
                        dtype={"CONTIG": str, "START": np.int64, "END": np.int64, "COUNT": np.float64})
    missing = [c for c in COORDS + ["COUNT"] if c not in table.columns]
    if missing:
        raise ValueError("%s: missing columns %s" % (path, missing))
    sample = _sample_from_header(header) or os.path.basename(path).split(".")[0]
    return sample, table[COORDS], table["COUNT"].to_numpy()


def find_count_files(directory: str, pattern: str = "*.tsv") -> List[str]:
    paths = sorted(p for p in glob.glob(os.path.join(directory, pattern)) if os.path.isfile(p))
    if not paths:
        raise FileNotFoundError("no files matching %s in %s" % (pattern, directory))
    return paths


def read_counts(paths: Sequence[str]) -> CountMatrix:
    """Read and align count files. All files must share identical intervals."""
    samples, rows = [], []
    ref = None
    for path in paths:
        sample, intervals, counts = read_count_file(path)
        if ref is None:
            ref = intervals.reset_index(drop=True)
        elif len(intervals) != len(ref) or not (
                np.array_equal(intervals["CONTIG"].to_numpy(), ref["CONTIG"].to_numpy())
                and np.array_equal(intervals["START"].to_numpy(), ref["START"].to_numpy())
                and np.array_equal(intervals["END"].to_numpy(), ref["END"].to_numpy())):
            raise ValueError("interval mismatch: %s does not have the same intervals as %s" % (path, paths[0]))
        if sample in samples:
            raise ValueError("duplicate sample name '%s' (%s)" % (sample, path))
        samples.append(sample)
        rows.append(counts)
    if ref is None:
        raise ValueError("no count files given")
    return CountMatrix(samples, ref, np.vstack(rows))
