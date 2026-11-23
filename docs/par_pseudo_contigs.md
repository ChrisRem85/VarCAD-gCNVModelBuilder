# PAR1/PAR2 as in-silico pseudo-autosomes (SHOX)

Design background: [gCNV_routing_audit_and_plan.md](gCNV_routing_audit_and_plan.md), section 1.3a. This page describes what is implemented in the COHORT workflow ([snakemake/gatk_gCNV.cohort.smk](../snakemake/gatk_gCNV.cohort.smk)) and how to use it.

## Purpose

GATK gCNV has no PAR option. On a GRCh38 reference with chrY PAR hard-masked, all PAR reads align to chrX PAR, so chrX PAR shows 2 copies in both 46,XX and 46,XY while the chrX contig ploidy is 1 in males. Left alone this gives a false DUP in every male and inflates the male chrX depth. Excluding PAR would make SHOX (PAR1) deletions uncallable.

Solution: chrX PAR1/PAR2 intervals are relabeled to two extra diploid contigs, `chrPAR1` and `chrPAR2`, before `FilterIntervals`. Everything downstream (ploidy, gCNV) treats them as autosomes. After postprocessing the calls are mapped back to chrX. Coordinates are never shifted: this is a pure contig rename.

## Config key `par_mode`

```bash
snakemake -s snakemake/gatk_gCNV.cohort.smk --config protocol=wgs.1k par_mode=pseudo_contig ...
```

| value | behaviour |
|-------|-----------|
| `none` (default) | Existing behaviour. No PAR rule is requested; all rule commands are the same as before. |
| `pseudo_contig` | PAR1/PAR2 become `chrPAR1`/`chrPAR2` (see below). For WGS/WES (SHOX). |
| `exclude` | Every chrX and chrY interval overlapping a PAR is dropped (implemented with `par_relabel.py --mode exclude`, not `-XL`, so counts, preprocessed and annotated intervals stay consistent). Dictionary, priors, scatter and postprocessing are unchanged. For panels without PAR targets. |

Optional: `par_exclude_bed=/path/to.bed` (any mode except `none`) additionally drops intervals overlapping a BED, for example the X-transposed region (XTR). The BED is 0-based half-open with original contig names, and the path must be inside a mounted directory (`CWD`, `SCRIPT_DIR` or `DB_DIR`). No XTR coordinates are hard-coded.

PAR coordinates used (GRCh38, 1-based closed): chrX PAR1 10001-2781479, PAR2 155701383-156030895; chrY PAR1 10001-2781479, PAR2 56887903-57217415.

## What `pseudo_contig` does

- Pseudo dictionary: the original `@SQ` list with `chrPAR1`, `chrPAR2` inserted immediately before chrX (right after chr22). `LN` is the chrX length. M5/UR of the originals are kept, the new contigs have no M5.
- chrX intervals fully inside PAR1/PAR2 are relabeled, intervals straddling a PAR boundary are dropped, all chrY intervals overlapping a chrY PAR are dropped.
- Rows are re-sorted in pseudo-dictionary order, the `@SQ` header is replaced, other header lines (`@HD`, `@RG`, ...) are kept.
- Priors: [assets/contig_ploidy_prior.par_pseudo.tsv](../assets/contig_ploidy_prior.par_pseudo.tsv) = the standard priors plus `chrPAR1`/`chrPAR2` `0.01 0.02 0.95 0.02`. The standard file is unchanged.
- Sharding: `chrPAR1` and `chrPAR2` go into the chr22 shard (`scatter_filtered_intervals`, same shard name, dictionary order).
- Postprocessing: `chrX`/`chrY` stay the only `--allosomal-contig`. `--sequence-dictionary <pseudo.dict>` replaces `--reference`.

## File flow

| rule | in | out |
|------|----|-----|
| `par_make_dict` | `assets/{protocol}/hg38.preprocessed.interval_list` (`@SQ`) | `gCNV/02a_par_pseudo_intervals/{protocol}/hg38.par_pseudo.dict` |
| `par_relabel_intervals` (`kind` = preprocessed, annotated) | `assets/{protocol}/hg38.{kind}.interval_list` | `gCNV/02a_par_*_intervals/{protocol}/hg38.{kind}.interval_list` |
| `par_relabel_read_counts` (per dataset) | `gCNV/03_read_counts/{protocol}/{dataset}.hg38.tsv` | `gCNV/03a_par_*_read_counts/{protocol}/{dataset}.hg38.tsv` |
| `filter_intervals` | relabeled counts + relabeled interval lists | `04_filtered_intervals` (unchanged path) |
| `determine_germline_contig_ploidy` | relabeled counts, `par_pseudo` priors (pseudo mode) | `05_contig_ploidy` |
| `scatter_filtered_intervals` | filtered list | `06_scattered_intervals`; chr22 shard also holds chrPAR1/2 |
| `run_GermlineCNVCaller` | relabeled counts, relabeled annotated intervals | `07_raw_cnv` |
| `postprocess_GermlineCNVCalls` | as before, `--sequence-dictionary` | `08_postprocessed_cnv` (pseudo contigs, raw GATK output) |
| `par_reconstruct` (per sample) | `08_postprocessed_cnv`, pseudo dict, sample `contig_ploidy.tsv` | `gCNV/08a_par_reconstructed/{protocol}/cohort.SAMPLE_{i}/` `genotyped_segments.vcf.gz`, `genotyped_intervals.vcf.gz`, `denoised_copy_ratios.tsv` |
| `par_ploidy_qc` (per sample) | sample `contig_ploidy.tsv` | `08a_par_reconstructed/.../par_ploidy_qc.tsv` |
| `finalize_genotyped_segments` | reconstructed segments (pseudo mode) | `gCNV/final/...genotyped_segments.vcf.gz` |

In `exclude` mode the directory tag is `par_exclude`, and only the relabel rules, `filter_intervals`, DGCP and gCNV are affected.

The scripts use only the Python 3 standard library (3.6 compatible) and run in the `broadinstitute/gatk` container via `python3`.

## Script usage

All scripts are in [scripts/par/](../scripts/par/) and have `--help`.

```bash
# pseudo dictionary from any file with an @SQ header
python3 scripts/par/par_make_dict.py --input hg38.preprocessed.interval_list --output hg38.par_pseudo.dict

# relabel counts TSV, annotated intervals TSV or interval_list (format is detected)
python3 scripts/par/par_relabel.py --input S1.hg38.tsv --output S1.hg38.par.tsv \
    [--mode pseudo_contig|exclude] [--exclude-bed xtr.bed]

# after PostprocessGermlineCNVCalls
python3 scripts/par/par_reconstruct.py --pseudo-dict hg38.par_pseudo.dict \
    --contig-ploidy contig_ploidy.tsv \
    --segments-vcf seg.vcf.gz   --out-segments-vcf seg.chrX.vcf \
    --intervals-vcf int.vcf.gz  --out-intervals-vcf int.chrX.vcf \
    --denoised-copy-ratios dcr.tsv --out-denoised-copy-ratios dcr.chrX.tsv \
    [--boundary-tolerance 0]
bcftools sort -O z -o seg.chrX.vcf.gz seg.chrX.vcf && bcftools index -t seg.chrX.vcf.gz

# ploidy consistency QC
python3 scripts/par/par_ploidy_qc.py --contig-ploidy contig_ploidy.tsv --output par_ploidy_qc.tsv
```

### Reconstruction details

- `chrPAR1`/`chrPAR2` become `chrX`; VCF `##contig` lines and TSV `@SQ` lines are replaced by the original dictionary (pseudo dictionary minus the pseudo contigs, or `--original-dict`). VCF `##contig` lines contain only `ID` and `length`.
- Records are sorted by the original dictionary (PAR1 at the start of chrX, PAR2 at the end). The VCFs are written as plain text, the rule sorts again with `bcftools sort -O z` and indexes with `bcftools index -t`.
- INFO tags (segments and intervals VCF): `PAR=PAR1|PAR2` on PAR records only; `REFCN=2` on PAR records, and on chrX-specific records the sample's chrX ploidy from `contig_ploidy.tsv`. chrY records get no tag.
- `PAR_BOUNDARY_ADJACENT` (flag, segments VCF only, non-reference segments only; non-reference = `ALT != '.'` and a genotype that is not all-0, as in `finalize_genotyped_segments`). Boundary positions come from the sample's intervals VCF, i.e. from the intervals that were actually modelled:
  - non-ref `chrPAR1` segment whose END is at or after the end of the last `chrPAR1` interval
  - non-ref chrX segment whose POS is at or before the start of the first chrX interval (first X-specific interval after PAR1)
  - non-ref `chrPAR2` segment whose POS is at or before the start of the first `chrPAR2` interval
  - non-ref chrX segment whose END is at or after the end of the last chrX interval (last X-specific interval before PAR2)

  `--boundary-tolerance` (bp, default 0) loosens each comparison. Segments are never merged; an analyst reports a flagged pair as one event. Denoised copy ratios are not rescaled to `REFCN`.

## QC output

`par_ploidy_qc.tsv` per sample, columns `sample contig observed expected PLOIDY_GQ status`. Expected ploidy(chrPAR1) = ploidy(chrX) + ploidy(chrY), same for chrPAR2.

```
sample  contig   observed  expected  PLOIDY_GQ  status
S1      chrPAR1  1         2         80.5       MISMATCH
```

`status` is `OK`, `MISMATCH` or `NA` (chrX/chrY/chrPARn missing in the ploidy calls). A mismatch is a finding (whole-PAR copy-number change, sample swap, sex-check problem), never a job failure. Examples: 45,X with X=1, Y=0, PAR=1 is OK; 46,XY with PAR=1 is a MISMATCH.

## Tests

```bash
python -m pytest tests/par -q
```

Synthetic files only; no GATK or Snakemake needed.

## To verify on the cluster

1. chrY PAR is hard-masked in the reference (0 = masked):
   `samtools faidx REF.fasta chrY:10001-2781479 | grep -v '>' | tr -d 'Nn\n' | wc -c`
2. `FilterIntervals`, `DetermineGermlineContigPloidy` and `GermlineCNVCaller` accept the relabeled counts, interval lists and annotated intervals (identical dictionaries, pseudo contigs in a sorted position, priors found for `chrPAR1`/`chrPAR2`).
3. `PostprocessGermlineCNVCalls` accepts `--sequence-dictionary <pseudo.dict>` without `--reference`, and produces correct VCF headers (`##contig` for chrPAR1/2, REF base handling, `END`, and `GT` of chrPAR records).
4. Behaviour for a 45,X sample: PAR ploidy 1 in the QC file, CN of PAR records, and how postprocessing labels a contig with ploidy 1.
5. The median male/female PAR1 depth ratio is about 1.0 (X-specific bins about 0.5).
6. `python3` in the `broadinstitute/gatk` image (the GATK conda env is Python 3.6) runs the scripts; `bcftools sort`/`bcftools index` exist there.
7. The pseudo-contig counts files are small enough for the `par_relabel_read_counts` memory (8 GB per sample, 16 GB for the interval lists) with your bin size.
8. Reconstructed VCFs index and load (`bcftools view`, IGV) and `finalize_genotyped_segments` still behaves with the extra INFO tags.

## Not covered here

- `build_gCNVModel.sh` (model build) is not changed. To build a PAR-aware model there it needs the same steps: relabel the read counts and interval lists, the `par_pseudo` priors, the chr22 shard including chrPAR1/2, and later the same postprocessing and reconstruction.
- CASE mode and the routing workflows (`route.smk`, `call.smk`) are not implemented in this repository yet; the same scripts apply.
- Unmasked chrY PAR (summing X + Y reads per bin) is not supported.
