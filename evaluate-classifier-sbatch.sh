#!/bin/bash
#SBATCH --job-name=evaluate-classifier
#SBATCH --time=48:00:00
#SBATCH --mem=100G

set -euo pipefail

log_status() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] [evaluate-classifier] $*"
}

usage() {
  cat <<'EOF'
Usage:
  sbatch evaluate-classifier-sbatch.sh \
    --classifier NEW_CLASSIFIER.qza \
    --old-taxonomy OLD_TAXONOMY.qza \
    --sequences REP_SEQS.qza \
    --table TABLE.qza \
    --output-dir OUTPUT_DIR \
    [--levels 2,3,4,5,6,7] \
    [--confidence 0.7]

Required arguments:
  --classifier    Path to FeatureData[Classifier] artifact (.qza)
  --old-taxonomy  Path to old FeatureData[Taxonomy] artifact (.qza)
  --sequences     Path to FeatureData[Sequence] artifact (.qza)
  --table         Path to FeatureTable[Frequency] artifact (.qza)
  --output-dir    Directory for evaluation outputs/report

Optional arguments:
  --levels        Comma-separated taxonomic levels (default: 2,3,4,5,6,7)
  --confidence    classify-sklearn confidence value (default: 0.7)
  -h, --help      Show this help message
EOF
}

CLASSIFIER=""
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

if [[ -z "$CLASSIFIER" || -z "$OLD_TAXONOMY" || -z "$SEQUENCES" || -z "$TABLE" || -z "$OUTPUT_DIR" ]]; then
  log_status "Missing required argument(s)."
  usage
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EVAL_SCRIPT="${SCRIPT_DIR}/evaluate_classifier.py"

if [[ ! -f "$EVAL_SCRIPT" ]]; then
  log_status "Could not find evaluation script: ${EVAL_SCRIPT}"
  exit 1
fi

log_status "Starting classifier evaluation"
log_status "Classifier: ${CLASSIFIER}"
log_status "Old taxonomy: ${OLD_TAXONOMY}"
log_status "Sequences: ${SEQUENCES}"
log_status "Table: ${TABLE}"
log_status "Output dir: ${OUTPUT_DIR}"
log_status "Levels: ${LEVELS}"
log_status "Confidence: ${CONFIDENCE}"

python "$EVAL_SCRIPT" \
  --classifier "$CLASSIFIER" \
  --old-taxonomy "$OLD_TAXONOMY" \
  --sequences "$SEQUENCES" \
  --table "$TABLE" \
  --output-dir "$OUTPUT_DIR" \
  --levels "$LEVELS" \
  --confidence "$CONFIDENCE"

log_status "Classifier evaluation complete"
