# Decouple insertion-level weight computation into its own rule

**Date:** 2026-08-18
**Status:** approved, ready to implement
**Scope:** pure refactor — numerically identical output

## Problem

Weight computation is spread across three places:

| Step | Location |
|---|---|
| `1 - clip(R²)` (no-replicates only) | `compute_r2_weights.py` |
| `-log10(clip(w))` (both branches) | `gene_level_depletion_analysis.py:135` |
| Gene-timepoint normalisation | `gene_level_depletion_analysis.py:194` |

`gene_level_depletion_analysis.py` therefore has to know which branch produced
its weights, and any change to the weighting algorithm means editing the
aggregation script. The goal is a single upstream rule that owns weights end to
end, so the aggregation script consumes a finished weight column and new
weighting schemes can be added without touching it.

## Current behaviour to preserve

All three quirks below are held up by `clip(upper=1-1e-6)` and are easy to lose
in a rewrite. Verified against `HD_DIT_HAP` (replicates) and `LD_haploid`
(no replicates).

1. **`padj == 1` does not become a zero weight.** In `HD_DIT_HAP`, 100% of
   in-gene `YES0` and 68% of `YES3` adjusted p-values are exactly 1 (BH
   ceiling). The upper clip maps all of them to the same `4.34e-7`, so those
   gene-timepoints collapse to an *arithmetic mean* rather than `0/0 = NaN`.
   Zero of 22,565 gene-timepoint groups have a zero weight sum, which is why the
   current code needs no uniform fallback.
2. **Negative R² takes the same path.** `LD_haploid` has 67,237 insertions with
   R² < 0 (fit worse than a flat line); all clip to the `1e-6` lower bound and
   receive the floor weight.
3. **NaN is caught twice.** 21 insertions with NaN R² produce 168 NaN cells in
   the weights file; `fillna(1)` then the upper clip give them the same
   `4.34e-7` floor.

`fillna(1)` must run **before** `clip`. Reversed, NaN bypasses the ceiling and
`-log10(1)` yields a 0 weight, so those cells stop voting entirely.

Normalisation is currently a numerical no-op — `np.average(x, weights=w)`
already computes `Σwx/Σw`. It is kept because it is the mount point for a
zero-sum fallback once gated schemes are added.

Normalisation happens *after* the LFC merge, and `stack()` drops NaN, so the
denominator only covers cells with an observed LFC. The new rule needs
`LFC.tsv` as an input to reproduce this.

## Signal availability per branch

| Branch | Available | Feasible schemes |
|---|---|---|
| Replicates (DESeq2) | `padj`, `lfcSE`, `baseMean`, `stat` | naive, and later A1/D2/D5 + gates |
| No replicates | `LFC`, `baseMean`, rule 15 `R²` | R²-based only |

The no-replicates branch has no p-values at all, so A1/D2/D5 are not
implementable there. This is a data constraint, not a design choice.

## Design

### New script

`workflow/scripts/depletion_scoring/compute_insertion_weights.py` replaces
`compute_r2_weights.py`, which is deleted.

```
load → in-gene filter → raw_weights(scheme) → mask(lfc.notna())
     → stack to long → join Systematic ID → normalise per gene-timepoint → write
```

```python
case Scheme.NAIVE:   # replicates
    weights = -np.log10(padj.fillna(1).clip(1e-6, 1-1e-6))
case Scheme.R2:      # no replicates
    conf = r2.clip(1e-6, 1-1e-6)
    weights = -np.log10((1 - conf).fillna(1).clip(1e-6, 1-1e-6))
```

`Scheme` is an enum so A1/D2/D5 can be added as new branches. `GateStyle` /
`GateSpec` are deliberately not introduced yet — there is no gate today.

### Output format

`16_gene_level_depletion_analysis/insertion_weights.tsv`, long format:
`Chr, Coordinate, Strand, Target, Timepoint, Systematic ID, Weight`.

Long rather than wide because `annotations.tsv` has 2,478 duplicate index rows
(one insertion annotated to two genes) and the normalising denominator is
gene-specific, so such an insertion must carry a different weight per gene. A
wide layout cannot express this. Cost: ~360–390k rows versus the current
90–150k.

This path replaces the existing `transformed_weights.tsv` — same directory, same
role — and is a *declared* rule output, so one entry drops off the anchor
workaround list in `packaging.smk`.

### Rule wiring

Both branches write the same output path; the branch decision moves up into the
weight rule, so `gene_level_depletion_analysis`'s `weights_path` `branch()` call
is deleted and it reads one fixed path.

```python
if config.get("use_DEseq2_for_biological_replicates", False):
    rule compute_insertion_weights:
        input:
            stats=".../14_.../padj.tsv",
            lfc=".../14_.../LFC.tsv",
            annotations=rules.concat_counts_and_annotations.output.annotations,
        params: scheme="naive"
else:
    rule compute_insertion_weights:
        input:
            stats=rules.insertion_level_curve_fitting.output.stats,
            lfc=".../14_.../LFC.tsv",
            annotations=rules.concat_counts_and_annotations.output.annotations,
        params: scheme="r2"
```

Neither `14_` nor `15_` works as the output directory: the no-replicates branch
depends on rule 15's output while the replicates branch depends only on rule 14,
so either choice misnumbers one branch.

`gene_level_depletion_analysis.py` loses the `-log10` transform,
`filter_in_gene_data`, and the weight half of `prepare_weighted_data`; it joins
the pre-computed long weights instead.

### Release contract

`release/insertion_level/transformed_weights.tsv` is replaced by
`insertion_level/insertion_weights.tsv`. The downstream streamlit app must be
updated to match.

## Verification

- Back up `16_.../LFC.tsv` for `HD_DIT_HAP` and `LD_haploid`, re-run, diff
  byte-for-byte. Must be empty. Both were reproduced during design with
  `max abs diff: 0.0`.
- Unit tests on both `raw_weights` branches covering `padj == 1`, negative R²,
  and NaN.

## Stale files (not touched by this work)

`16_.../gene_level_statistics.tsv` (2026-07-18) and `LFC.tsv` (2026-08-17) are
written from the same variable at lines 351/353 yet have opposite signs — the
former is a leftover from a run predating the sign fix in `1d36ee2`. `LFC2.tsv`
in the same directory is referenced by no rule. Use `LFC.tsv` as the
verification baseline.
