#!/bin/bash
#SBATCH --job-name=evaluate-classifier
#SBATCH --time=48:00:00
#SBATCH --mem=100G
#SBATCH --partition rhel10

set -euo pipefail

log_status() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] [evaluate-classifier] $*"
}

usage() {
  cat <<'EOF'
Usage:
  sbatch evaluate-classifier-sbatch.sh \
    (--classifier NEW_CLASSIFIER.qza --sequences REP_SEQS.qza | --new-taxonomy NEW_TAXONOMY.qza) \
    --old-taxonomy OLD_TAXONOMY.qza \
    --table TABLE.qza \
    --output-dir OUTPUT_DIR \
    [--levels 2,3,4,5,6,7] \
    [--confidence 0.7]

Required arguments (mutually exclusive):
  --classifier    Path to FeatureData[Classifier] artifact (.qza); requires --sequences
  --new-taxonomy  Path to a pre-computed FeatureData[Taxonomy] artifact (.qza);
                  skips classify-sklearn

Other required arguments:
  --old-taxonomy  Path to old FeatureData[Taxonomy] artifact (.qza)
  --table         Path to FeatureTable[Frequency] artifact (.qza)
  --output-dir    Directory for evaluation outputs/report

Optional arguments:
  --sequences     Path to FeatureData[Sequence] artifact (.qza); required with --classifier
  --levels        Comma-separated taxonomic levels (default: 2,3,4,5,6,7)
  --confidence    classify-sklearn confidence value (default: 0.7)
  -h, --help      Show this help message
EOF
}

CLASSIFIER=""
NEW_TAXONOMY=""
OLD_TAXONOMY=""
SEQUENCES=""
TABLE=""
OUTPUT_DIR=""
LEVELS="2,3,4,5,6,7"
CONFIDENCE="0.7"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --classifier)
      CLASSIFIER="$2"
      shift 2
      ;;
    --new-taxonomy)
      NEW_TAXONOMY="$2"
      shift 2
      ;;
    --old-taxonomy)
      OLD_TAXONOMY="$2"
      shift 2
      ;;
    --sequences)
      SEQUENCES="$2"
      shift 2
      ;;
    --table)
      TABLE="$2"
      shift 2
      ;;
    --output-dir)
      OUTPUT_DIR="$2"
      shift 2
      ;;
    --levels)
      LEVELS="$2"
      shift 2
      ;;
    --confidence)
      CONFIDENCE="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      log_status "Unknown argument: $1"
      usage
      exit 1
      ;;
  esac
done

if [[ -z "$OLD_TAXONOMY" || -z "$TABLE" || -z "$OUTPUT_DIR" ]]; then
  log_status "Missing required argument(s)."
  usage
  exit 1
fi

if [[ -z "$CLASSIFIER" && -z "$NEW_TAXONOMY" ]]; then
  log_status "Provide either --classifier or --new-taxonomy."
  usage
  exit 1
fi

if [[ -n "$CLASSIFIER" && -n "$NEW_TAXONOMY" ]]; then
  log_status "--classifier and --new-taxonomy are mutually exclusive."
  usage
  exit 1
fi

if [[ -n "$CLASSIFIER" && -z "$SEQUENCES" ]]; then
  log_status "--sequences is required when --classifier is used."
  usage
  exit 1
fi

# Locate evaluate_classifier.py next to this script. Under Slurm the batch
# script is executed from a spool copy, so BASH_SOURCE no longer points at the
# repo; recover the submitted script's path from scontrol in that case.
if [[ -n "${SLURM_JOB_ID:-}" ]] && command -v scontrol >/dev/null 2>&1; then
  SBATCH_SCRIPT="$(scontrol show job "$SLURM_JOB_ID" 2>/dev/null | sed -n 's/^ *Command=//p' | head -1)"
  SBATCH_SCRIPT="${SBATCH_SCRIPT%% *}"
  SCRIPT_DIR="$(cd "$(dirname "$SBATCH_SCRIPT")" && pwd)"
else
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
EVAL_SCRIPT="${SCRIPT_DIR}/evaluate_classifier.py"

if [[ ! -f "$EVAL_SCRIPT" ]]; then
  log_status "Could not find evaluation script: ${EVAL_SCRIPT}"
  exit 1
fi

log_status "Starting classifier evaluation"
log_status "Classifier: ${CLASSIFIER:-<not provided>}"
log_status "New taxonomy: ${NEW_TAXONOMY:-<not provided>}"
log_status "Old taxonomy: ${OLD_TAXONOMY}"
log_status "Sequences: ${SEQUENCES:-<not provided>}"
log_status "Table: ${TABLE}"
log_status "Output dir: ${OUTPUT_DIR}"
log_status "Levels: ${LEVELS}"
log_status "Confidence: ${CONFIDENCE}"

TAXONOMY_ARG=""
if [[ -n "$CLASSIFIER" ]]; then
  TAXONOMY_ARG="--classifier ${CLASSIFIER} --sequences ${SEQUENCES}"
else
  TAXONOMY_ARG="--new-taxonomy ${NEW_TAXONOMY}"
fi

# shellcheck disable=SC2086
python "$EVAL_SCRIPT" \
  $TAXONOMY_ARG \
  --old-taxonomy "$OLD_TAXONOMY" \
  --table "$TABLE" \
  --output-dir "$OUTPUT_DIR" \
  --levels "$LEVELS" \
  --confidence "$CONFIDENCE"

log_status "Classifier evaluation complete"
