#!/bin/bash
#SBATCH --job-name=evaluate-classifier
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --partition rhel10

# Compare two taxonomy assignments ("old" and "new") over the same sequences
# and feature table.
#
# Both taxonomies must already exist. Classification is deliberately out of
# scope: run classify-sklearn separately and pass the resulting artifact in.
# That keeps this script cheap to re-run, since the expensive, memory-hungry
# part of the workflow is no longer coupled to the comparison.
#
# Every analysis step below is a stock QIIME 2 command; this script only wires
# them together and bundles the resulting .qzv files into a single report.
# There is no custom Python analysis code.
#
# Requires an activated QIIME 2 2026.7 environment (e.g. rachis-qiime2-2026.7).

set -euo pipefail

log_status() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] [evaluate-classifier] $*"
}

warn() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] [evaluate-classifier] WARNING: $*" >&2
}

usage() {
  cat <<'EOF'
Usage:
  sbatch evaluate-classifier-sbatch.sh \
    --old-taxonomy OLD_TAXONOMY.qza \
    --new-taxonomy NEW_TAXONOMY.qza \
    --sequences REP_SEQS.qza \
    --table TABLE.qza \
    --output-dir OUTPUT_DIR \
    [--old-label old] [--new-label new] \
    [--levels 2,3,4,5,6,7]

Required arguments:
  --old-taxonomy    Path to old FeatureData[Taxonomy] artifact (.qza)
  --new-taxonomy    Path to new FeatureData[Taxonomy] artifact (.qza)
  --sequences       Path to the classified FeatureData[Sequence] artifact (.qza)
  --table           Path to FeatureTable[Frequency] artifact (.qza)
  --output-dir      Directory for evaluation outputs/report

Optional arguments:
  --old-label       Label for the old taxonomy in visualizations
                    (default: the --old-taxonomy file name, minus .qza)
  --new-label       Label for the new taxonomy in visualizations
                    (default: the --new-taxonomy file name, minus .qza)
  --levels          Comma-separated taxonomic levels to collapse/compare (default: 2,3,4,5,6,7)
  -h, --help        Show this help message

To produce a new taxonomy for comparison, run classification first:
  qiime feature-classifier classify-sklearn \
    --i-classifier CLASSIFIER.qza \
    --i-reads REP_SEQS.qza \
    --p-confidence 0.7 \
    --o-classification NEW_TAXONOMY.qza
EOF
}

NEW_TAXONOMY=""
OLD_TAXONOMY=""
SEQUENCES=""
TABLE=""
OUTPUT_DIR=""
OLD_LABEL=""
NEW_LABEL=""
LEVELS="2,3,4,5,6,7"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --new-taxonomy)    NEW_TAXONOMY="$2"; shift 2 ;;
    --old-taxonomy)    OLD_TAXONOMY="$2"; shift 2 ;;
    --sequences)       SEQUENCES="$2"; shift 2 ;;
    --table)           TABLE="$2"; shift 2 ;;
    --output-dir)      OUTPUT_DIR="$2"; shift 2 ;;
    --old-label)       OLD_LABEL="$2"; shift 2 ;;
    --new-label)       NEW_LABEL="$2"; shift 2 ;;
    --levels)          LEVELS="$2"; shift 2 ;;
    -h|--help)         usage; exit 0 ;;
    *)                 log_status "Unknown argument: $1"; usage; exit 1 ;;
  esac
done

if [[ -z "$OLD_TAXONOMY" || -z "$NEW_TAXONOMY" || -z "$SEQUENCES" || -z "$TABLE" || -z "$OUTPUT_DIR" ]]; then
  log_status "Missing required argument(s)."
  usage
  exit 1
fi

if ! command -v qiime >/dev/null 2>&1; then
  log_status "The 'qiime' command was not found; activate a QIIME 2 2026.7 environment first."
  exit 1
fi

# Label the two taxonomies by their file names unless the caller says otherwise.
# These labels appear in the rescript summary plots, the concordance plot, and
# both axes of every Mantel result, so defaulting them to the artifact names
# keeps a report self-describing: months later it is still obvious which two
# taxonomies produced it.
if [[ -z "$OLD_LABEL" ]]; then
  OLD_LABEL="$(basename "$OLD_TAXONOMY" .qza)"
fi
if [[ -z "$NEW_LABEL" ]]; then
  NEW_LABEL="$(basename "$NEW_TAXONOMY" .qza)"
fi

# Parse --levels into an array, rejecting anything that is not a usable depth
# before we spend time on the earlier steps.
IFS=',' read -r -a LEVEL_ARRAY <<< "$LEVELS"
for level in "${LEVEL_ARRAY[@]}"; do
  if ! [[ "$level" =~ ^[0-9]+$ ]] || (( level < 1 )); then
    log_status "Invalid taxonomic level '${level}'; levels must be positive integers."
    exit 1
  fi
done

ARTIFACT_DIR="${OUTPUT_DIR}/artifacts"
VIZ_DIR="${OUTPUT_DIR}/visualizations"
mkdir -p "$ARTIFACT_DIR" "$VIZ_DIR"

# Visualizations accumulate here and are bundled into one report at the end.
# Steps that fail are skipped rather than aborting the run, so a single
# unsupported comparison does not cost you the whole job.
REPORT_VIZ=()

log_status "Starting classifier evaluation"
log_status "Old taxonomy:  ${OLD_TAXONOMY} (labelled '${OLD_LABEL}')"
log_status "New taxonomy:  ${NEW_TAXONOMY} (labelled '${NEW_LABEL}')"
log_status "Sequences:     ${SEQUENCES}"
log_status "Table:         ${TABLE}"
log_status "Output dir:    ${OUTPUT_DIR}"
log_status "Levels:        ${LEVELS}"

# ---------------------------------------------------------------------------
# Step 1 - Summary statistics for each taxonomy on its own.
#
# Reports unique labels, taxonomic entropy, and the number of (un)classified
# features at each level, plotted for both taxonomies side by side. This is
# descriptive, not a comparison: it answers "how deep and how diverse is each
# assignment" before we ask how much they agree.
# ---------------------------------------------------------------------------
log_status "Step 1: rescript evaluate-taxonomy (per-taxonomy summary stats)"
if qiime rescript evaluate-taxonomy \
     --i-taxonomies "$OLD_TAXONOMY" \
     --i-taxonomies "$NEW_TAXONOMY" \
     --p-labels "$OLD_LABEL" \
     --p-labels "$NEW_LABEL" \
     --o-taxonomy-stats "${VIZ_DIR}/rescript-evaluate-taxonomy.qzv"; then
  REPORT_VIZ+=("${VIZ_DIR}/rescript-evaluate-taxonomy.qzv")
else
  warn "rescript evaluate-taxonomy failed; continuing."
fi

# ---------------------------------------------------------------------------
# Step 2 - Per-level concordance between the two taxonomies.
#
# Treats the old taxonomy as "expected" and the new one as "observed" and plots
# precision, recall, and F-measure at each level. Note this measures agreement,
# not correctness: the old assignment is a reference point, not ground truth.
# ---------------------------------------------------------------------------
log_status "Step 2: rescript evaluate-classifications (old-vs-new concordance)"
if qiime rescript evaluate-classifications \
     --i-expected-taxonomies "$OLD_TAXONOMY" \
     --i-observed-taxonomies "$NEW_TAXONOMY" \
     --p-labels "${OLD_LABEL}-vs-${NEW_LABEL}" \
     --o-evaluation "${VIZ_DIR}/rescript-evaluate-classifications.qzv"; then
  REPORT_VIZ+=("${VIZ_DIR}/rescript-evaluate-classifications.qzv")
else
  warn "rescript evaluate-classifications failed; continuing."
fi

# ---------------------------------------------------------------------------
# Step 3 - Feature-by-feature taxonomy comparison table.
#
# Puts both taxonomy strings side by side on one sortable row per sequence,
# alongside the least common ancestor of the two. The LCA truncates at the
# level where the assignments stop agreeing, so reading across a row shows both
# what each classifier called the feature and how deep the agreement ran.
#
# tabulate-seqs takes its taxonomies as a Collection, whose keys become the
# column headers. The CLI builds a Collection from a directory, keyed by file
# name, with an optional .order file fixing column order. Note the keys are
# restricted to [A-Za-z0-9+-._~] -- no spaces -- so the agreement column is
# "Classifier-agreement" rather than "Classifier agreement".
#
# merge-method 'intersect' restricts the table to features that appear in the
# feature table as well as the sequences and taxonomies. Reference taxonomies
# routinely cover far more features than a study's table does, and without this
# the visualization would be dominated by rows carrying no abundance data. Use
# 'union' instead to tabulate every classified sequence.
# ---------------------------------------------------------------------------
log_status "Step 3: feature-by-feature taxonomy comparison (tabulate-seqs)"
LCA_TAXONOMY="${ARTIFACT_DIR}/classifier-agreement-lca.qza"
FEATURE_FREQS="${ARTIFACT_DIR}/feature-frequencies.qza"
TAXA_COLLECTION="${ARTIFACT_DIR}/taxonomy-collection"

if qiime rescript merge-taxa \
     --i-data "$OLD_TAXONOMY" \
     --i-data "$NEW_TAXONOMY" \
     --p-mode lca \
     --o-merged-data "$LCA_TAXONOMY" \
  && qiime feature-table tabulate-feature-frequencies \
     --i-table "$TABLE" \
     --o-feature-frequencies "$FEATURE_FREQS"; then

  rm -rf "$TAXA_COLLECTION"
  mkdir -p "$TAXA_COLLECTION"
  cp "$LCA_TAXONOMY" "${TAXA_COLLECTION}/Classifier-agreement.qza"
  cp "$OLD_TAXONOMY" "${TAXA_COLLECTION}/${OLD_LABEL}.qza"
  cp "$NEW_TAXONOMY" "${TAXA_COLLECTION}/${NEW_LABEL}.qza"
  printf '%s\n%s\n%s\n' "Classifier-agreement" "$OLD_LABEL" "$NEW_LABEL" \
    > "${TAXA_COLLECTION}/.order"

  if qiime feature-table tabulate-seqs \
       --i-data "$SEQUENCES" \
       --i-taxonomy "$TAXA_COLLECTION" \
       --m-metadata-file "$FEATURE_FREQS" \
       --p-merge-method intersect \
       --o-visualization "${VIZ_DIR}/taxonomy-comparison.qzv"; then
    REPORT_VIZ+=("${VIZ_DIR}/taxonomy-comparison.qzv")
  else
    warn "feature-table tabulate-seqs failed; continuing."
  fi
else
  warn "Could not build the LCA taxonomy or feature frequencies; skipping the comparison table."
fi

# ---------------------------------------------------------------------------
# Step 4 - Does the choice of taxonomy change the community structure?
#
# For each requested level, collapse the same feature table twice (once per
# taxonomy), compute Bray-Curtis distances between samples from each collapsed
# table, and Mantel-test the two distance matrices against each other.
#
# This is the question that matters for downstream ecology: even where the two
# taxonomies disagree on labels, a high Mantel correlation means the
# sample-to-sample relationships you would infer are effectively unchanged.
#
# Correlation does not vary with level in any fixed direction. A disagreement
# only moves abundance between bins at levels where the two taxonomies group
# features differently, so which levels suffer depends on where the classifiers
# diverge. Read the levels as independent comparisons rather than a trend.
# ---------------------------------------------------------------------------
log_status "Step 4: per-level Bray-Curtis / Mantel comparison"
for level in "${LEVEL_ARRAY[@]}"; do
  log_status "  Level ${level}: collapsing, computing Bray-Curtis, running Mantel test"

  old_collapsed="${ARTIFACT_DIR}/collapsed-old-level-${level}.qza"
  new_collapsed="${ARTIFACT_DIR}/collapsed-new-level-${level}.qza"
  old_dm="${ARTIFACT_DIR}/braycurtis-old-level-${level}.qza"
  new_dm="${ARTIFACT_DIR}/braycurtis-new-level-${level}.qza"
  mantel_viz="${VIZ_DIR}/mantel-level-${level}.qzv"

  if qiime taxa collapse \
       --i-table "$TABLE" \
       --i-taxonomy "$OLD_TAXONOMY" \
       --p-level "$level" \
       --o-collapsed-table "$old_collapsed" \
    && qiime taxa collapse \
       --i-table "$TABLE" \
       --i-taxonomy "$NEW_TAXONOMY" \
       --p-level "$level" \
       --o-collapsed-table "$new_collapsed" \
    && qiime diversity beta \
       --i-table "$old_collapsed" \
       --p-metric braycurtis \
       --o-distance-matrix "$old_dm" \
    && qiime diversity beta \
       --i-table "$new_collapsed" \
       --p-metric braycurtis \
       --o-distance-matrix "$new_dm" \
    && qiime diversity mantel \
       --i-dm1 "$old_dm" \
       --i-dm2 "$new_dm" \
       --p-method spearman \
       --p-label1 "${OLD_LABEL} (level ${level})" \
       --p-label2 "${NEW_LABEL} (level ${level})" \
       --p-intersect-ids \
       --o-visualization "$mantel_viz"; then
    REPORT_VIZ+=("$mantel_viz")
  else
    warn "Level ${level} comparison failed (the taxonomy may not reach this depth); continuing."
  fi
done

# ---------------------------------------------------------------------------
# Step 5 - Bundle every visualization into a single browsable report.
#
# make-report writes to <report-path>.qzv, so the result is one file that can
# be opened in QIIME 2 View or shared as-is.
# ---------------------------------------------------------------------------
if (( ${#REPORT_VIZ[@]} == 0 )); then
  log_status "No visualizations were produced; skipping report generation."
  exit 1
fi

log_status "Step 5: bundling ${#REPORT_VIZ[@]} visualizations with tools make-report"
qiime tools make-report \
  --report-path "${OUTPUT_DIR}/evaluation-report" \
  "${REPORT_VIZ[@]}"

log_status "Classifier evaluation complete"
log_status "Report:         ${OUTPUT_DIR}/evaluation-report.qzv"
log_status "Visualizations: ${VIZ_DIR}"
log_status "Intermediates:  ${ARTIFACT_DIR}"
