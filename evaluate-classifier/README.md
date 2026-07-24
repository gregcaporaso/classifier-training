# evaluate-classifier

> **Status: experimental.**
> This is an internal, exploratory utility, not a supported or released tool.
> Interfaces, outputs, and file layout may change or be removed without notice.
> Do not depend on it in automated pipelines.

## What this is

A helper for sanity-checking a newly trained QIIME 2 taxonomic classifier by comparing its taxonomy assignments against an existing ("old") assignment on the *same* sequences and feature table.

It answers the question "how does the new classifier's output differ from what I had before?" and produces:

- **Feature-level taxonomy diff** — which features are assigned a different taxonomy, and the first taxonomic level at which the two assignments diverge.
- **Per-level collapsed-table comparison** — for each requested level, the counts of shared / old-only / new-only taxa, total frequencies, and the Spearman correlation of shared-taxon frequencies across samples.
- **Per-level frequency scatter plots** — log1p old vs. new counts for each (sample, taxon) pair among shared taxa.
- **A combined report** — the above visualizations bundled via `qiime tools make-report`.

The new taxonomy can be supplied directly (`--new-taxonomy`) or computed on the fly from a classifier (`--classifier` + `--sequences`, which runs `classify-sklearn`).

## What this is *not*

- **Not a ground-truth accuracy benchmark.** Every metric here is old-vs-new *concordance*. It tells you what changed and how well the two agree, not which classifier is objectively more correct. For accuracy against known answers, use a cross-validated benchmarking approach (e.g. `q2-feature-classifier` evaluation actions).
- **Not a QIIME 2 plugin.** `_actions.py` holds functions written to be *registration-ready* (each carries a registration sketch in its docstring), but they are currently imported and called directly by `evaluate_classifier.py`. This is a staging ground toward a possible future plugin, not a registered one.
- **Not a stable API.** The `_write_visualization` helper in `evaluate_classifier.py` hand-packages a `.qzv` and is an explicit stopgap; it goes away once `compare_collapsed_tables` becomes a real plugin visualizer.

## Files

- `evaluate_classifier.py` — the Click CLI entry point.
- `_actions.py` — the comparison functions (`compare_taxonomy`, `compare_collapsed_tables`, `build_scatter_metadata`), written for eventual plugin registration.
- `evaluate-classifier-sbatch.sh` — Slurm wrapper that invokes the CLI with the same arguments.

## Running it

This is memory-bound on the classifier: loading a full-length GTDB naive-Bayes model unpickles to several GiB, and `classify-sklearn` needs several times that at peak.
It will be killed by the OS (SIGKILL / out-of-memory) on a typical laptop; run it on an HPC node.
The provided Slurm wrapper requests `--mem=100G` for this reason.

```bash
sbatch evaluate-classifier-sbatch.sh \
  --classifier   NEW_CLASSIFIER.qza \
  --sequences    REP_SEQS.qza \
  --old-taxonomy OLD_TAXONOMY.qza \
  --table        TABLE.qza \
  --output-dir   OUTPUT_DIR \
  --levels 2,3,4,5,6,7 \
  --confidence 0.7
```

Pass `--new-taxonomy PRECOMPUTED_TAXONOMY.qza` in place of `--classifier`/`--sequences` to skip classification and compare an already-computed taxonomy.
