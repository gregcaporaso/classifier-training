# evaluate-classifier

> **Status: experimental.**
> This is an internal, exploratory utility, not a supported or released tool.
> Interfaces, outputs, and file layout may change or be removed without notice.
> Do not depend on it in automated pipelines.

## What this is

A helper for sanity-checking a newly trained QIIME 2 taxonomic classifier by comparing its taxonomy assignments against an existing ("old") assignment on the *same* sequences and feature table.

It answers the question "how does the new classifier's output differ from what I had before?"

The whole thing is a single Slurm batch script that calls stock QIIME 2 commands and bundles the resulting `.qzv` files into one report.
There is no custom analysis code: every number in the report is produced by a released QIIME 2 action.

Classification is out of scope.
Both taxonomies must already exist as `FeatureData[Taxonomy]` artifacts; run `classify-sklearn` yourself and pass the result in.
Keeping the expensive step separate means the comparison is cheap to re-run — you can change the levels or the labels and re-report in minutes without re-classifying.

## What it does

1. **Per-taxonomy summary stats** — `rescript evaluate-taxonomy` reports unique labels, taxonomic entropy, and the number of (un)classified features at each level, for both taxonomies side by side.
   This is descriptive rather than comparative: it shows how deep and how diverse each assignment is before asking how much the two agree.
2. **Per-level concordance** — `rescript evaluate-classifications` treats the old taxonomy as "expected" and the new one as "observed", plotting precision, recall, and F-measure at each level.
3. **Feature-by-feature comparison table** — `feature-table tabulate-seqs` puts both taxonomy strings side by side on one sortable row per sequence, together with the least common ancestor of the two (from `rescript merge-taxa --p-mode lca`) and the per-feature frequencies (from `feature-table tabulate-feature-frequencies`).
   The LCA truncates at the level where the two assignments stop agreeing, so reading across a row shows both what each classifier called the feature and how deep the agreement ran.
4. **Per-level community-structure comparison** — for each requested level, the feature table is collapsed twice (once per taxonomy) with `taxa collapse`, Bray-Curtis distances between samples are computed from each collapsed table with `diversity beta`, and the two distance matrices are compared with `diversity mantel`.
   This is the question that matters downstream: even where the labels disagree, a high Mantel correlation means the sample-to-sample relationships you would infer are effectively unchanged.
   Correlation does not vary with level in any fixed direction.
   A disagreement only moves abundance between bins at levels where the two taxonomies actually group features differently, so which levels suffer depends on where the classifiers diverge.
   Read the levels as independent comparisons rather than as a trend.
5. **Report** — `tools make-report` bundles every visualization into a single `evaluation-report.qzv`.

## What this is *not*

- **Not a ground-truth accuracy benchmark.**
  Every metric here is old-vs-new *concordance*.
  It tells you what changed and how well the two agree, not which classifier is objectively more correct.
  For accuracy against known answers, use a cross-validated benchmarking approach.
- **Not able to distinguish a rename from a reclassification.**
  When a reference database renames a taxon between releases, every affected feature gets a different taxonomy string even though nothing about its placement changed.
  Neither the string comparison in step 3 nor the Mantel tests in step 4 can tell that apart from a genuine reassignment, because both relabel bins identically.
  See [Future directions](#future-directions).
- **Not a QIIME 2 plugin.**
  It is a shell script that shells out to `qiime`.

## Files

- `evaluate-classifier-sbatch.sh` — the entire tool: Slurm wrapper, argument parsing, and the QIIME 2 pipeline.

## Running it

Requires an activated QIIME 2 environment (e.g. `rachis-qiime2-2026.7`).
This was tested with 2026.7, so that is presumed to be the lowest version that will work.

First produce the new taxonomy, if you do not already have it:

```bash
qiime feature-classifier classify-sklearn \
  --i-classifier NEW_CLASSIFIER.qza \
  --i-reads      REP_SEQS.qza \
  --p-confidence 0.7 \
  --o-classification NEW_TAXONOMY.qza
```

That step is the memory-hungry one: a full-length GTDB naive-Bayes model unpickles to several GiB and `classify-sklearn` peaks at a multiple of that, so it needs an HPC node with a large memory allocation.

Then run the comparison:

```bash
sbatch evaluate-classifier-sbatch.sh \
  --old-taxonomy 2026.4.0-r226-taxonomy.qza \
  --new-taxonomy 2026.7.0-r232-taxonomy.qza \
  --sequences    REP_SEQS.qza \
  --table        TABLE.qza \
  --output-dir   OUTPUT_DIR \
  --levels 2,3,4,5,6,7
```

`--sequences` is the artifact that was classified, and is what step 3 tabulates one row per.

`--old-label` and `--new-label` name the two taxonomies in the visualizations.
They default to the artifact file names with `.qza` stripped, so the example above labels its plots `2026.4.0-r226-taxonomy` and `2026.7.0-r226-taxonomy` without you having to say so.
The labels appear in the RESCRIPt summary plots, the concordance plot, and on both axes of every Mantel scatter plot, which keeps a report self-describing long after the run.

### Slurm is optional

Nothing in the script calls Slurm.
The `#SBATCH` lines are comments that `bash` ignores, so it runs locally just as well:

```bash
bash evaluate-classifier-sbatch.sh --old-taxonomy ... --new-taxonomy ... --sequences ... --table ... --output-dir ...
```

Local runs are practical now that classification is a separate step: the comparison itself is not memory-hungry, and every test of this script was run on a laptop.
Use `--levels 2` for a fast smoke test before committing to all six levels.

When you do submit it to Slurm, the script requests `--mem=64G` and `--time=12:00:00`.
These are a starting point rather than a measured requirement — `diversity beta` scales with the square of the sample count, so raise them for a large study.

## Outputs

- `OUTPUT_DIR/evaluation-report.qzv` — the combined report; open in QIIME 2 View.
- `OUTPUT_DIR/visualizations/` — the individual `.qzv` files that went into the report.
- `OUTPUT_DIR/artifacts/` — intermediate `.qza` files: the LCA "classifier agreement" taxonomy, per-feature frequencies, the taxonomy collection directory handed to `tabulate-seqs`, and the collapsed tables and distance matrices for each level.

## Notes

Individual steps are non-fatal.
A step that fails logs a warning and is left out of the report rather than aborting the run, so one unsupported comparison does not cost you a long job.
If a requested level is deeper than the taxonomy reaches, that level is skipped this way.

Expect a fixed overhead of roughly ten minutes independent of data size.
The script makes on the order of forty separate `qiime` invocations, each paying the QIIME 2 CLI startup cost.
That is the trade-off for using only released QIIME 2 actions rather than driving the framework from a single Python process.

For scale, a comparison of 1105 samples by 41 328 features, with 162 331 classified sequences, took 13.5 minutes end to end on a laptop and peaked at about 2 GB of memory.
Two thirds of that was step 4.
Bray-Curtis memory scales with the square of the *sample* count rather than the feature count, which is why a study with many features is cheaper than it looks.

## Future directions

### Distinguish reference renames from genuine reclassifications

The most useful thing this tool cannot currently do is tell you *why* a taxonomy string changed.

When comparing two releases of a reference database, a large share of the differing assignments are often pure nomenclature: the reference renamed a taxon, so every feature beneath it gets a new string while its phylogenetic placement is untouched.
A genuine reclassification — the classifier putting a sequence somewhere materially different — looks identical from where this script stands.
Both produce a differing string in step 3, and both relabel bins in step 4, so a rename depresses the Mantel correlation exactly as a reassignment would.
This matters most at shallow levels, where a single renamed phylum moves a large block of abundance between bins.

GTDB deliberately does not publish a release-to-release changelog that would resolve this.
The GTDB team's recommendation is a **genome-centric** comparison instead, on the grounds that taxa split, merge, and get reclassified in ways a rename list cannot represent.
The approach is to take the per-release taxonomy files (`bac120_taxonomy.tsv`, `ar53_taxonomy.tsv`) from each release and join them on **genome accession**.
Where the same accession carries a different lineage, the reference genuinely moved it; where the lineage differs only in the spelling of a rank, it was renamed.

Folding that into this report is not a matter of adding a metadata file.
`tabulate-seqs --m-metadata-file` needs metadata indexed by **feature ID**, whereas every GTDB artifact is keyed by taxon name or genome accession.
Bridging feature → taxon → release-diff is a three-way join, which is exactly the custom comparison logic this script exists to avoid.
It is better treated as a separate analysis whose output could then be passed in as a pre-computed, feature-indexed metadata file.

Relevant references:

- [GTDB forum: taxonomy name changes between releases](https://forum.gtdb.ecogenomic.org/t/taxonomy-name-changes-mapping-between-r220-and-r226/794) — the GTDB team's genome-centric recommendation.
- [GTDB release data](https://data.gtdb.ecogenomic.org/releases/) — per-release `*_taxonomy.tsv` files, plus `RELEASE_NOTES.txt` and `synonyms.*.tsv` in `auxillary_files/`.
- [GTDB Taxon History](https://gtdb.ecogenomic.org/taxon-history) — interactive per-taxon history across releases; useful for spot checks, not for bulk annotation.

### Other ideas

- **Report confidence shifts.**
  Nothing here looks at the `Confidence` column, so two taxonomies that assign identical strings with very different confidence appear perfectly concordant.
