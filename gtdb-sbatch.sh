#!/bin/bash
#SBATCH --job-name=train-gtdb
#SBATCH --time=48:00:00
#SBATCH --mem=100G

set -e

log_status() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] [train-gtdb] $*"
}

GTDB_VERSION="226.0"
GTDB_DOMAIN="Both"
GTDB_DB_TYPE="SpeciesReps"
GTDB_URL_TYPE="Primary"
GTDB_LABEL="gtdb_r${GTDB_VERSION}_${GTDB_DOMAIN}_${GTDB_DB_TYPE}"
GTDB_SEQS_ARTIFACT="${GTDB_LABEL}_seqs.qza"
GTDB_TAXONOMY_ARTIFACT="${GTDB_LABEL}_taxonomy.qza"
GTDB_CLASSIFIER_ARTIFACT="${GTDB_LABEL}_classifier.qza"

log_status "Pipeline starting for ${GTDB_LABEL}"
log_status "Checking whether classifier artifact already exists"

if [[ -f "${GTDB_CLASSIFIER_ARTIFACT}" ]]; then
  log_status "Classifier already exists: ${GTDB_CLASSIFIER_ARTIFACT}"
  log_status "Bailing out without downloading GTDB data or retraining"
  exit 0
fi

log_status "Classifier not found; checking for existing GTDB reference artifacts"

if [[ -f "${GTDB_SEQS_ARTIFACT}" && -f "${GTDB_TAXONOMY_ARTIFACT}" ]]; then
  log_status "Using existing GTDB artifacts: ${GTDB_SEQS_ARTIFACT} and ${GTDB_TAXONOMY_ARTIFACT}"
else
  log_status "Downloading GTDB reference data with RESCRIPt"
  qiime rescript get-gtdb-data \
    --p-version "${GTDB_VERSION}" \
    --p-domain "${GTDB_DOMAIN}" \
    --p-db-type "${GTDB_DB_TYPE}" \
    --p-url-type "${GTDB_URL_TYPE}" \
    --o-gtdb-sequences "${GTDB_SEQS_ARTIFACT}" \
    --o-gtdb-taxonomy "${GTDB_TAXONOMY_ARTIFACT}"
  log_status "GTDB reference artifacts created"
fi

log_status "Training naive Bayes classifier"
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads "${GTDB_SEQS_ARTIFACT}" \
  --i-reference-taxonomy "${GTDB_TAXONOMY_ARTIFACT}" \
  --o-classifier "${GTDB_CLASSIFIER_ARTIFACT}"

log_status "Classifier training complete: ${GTDB_CLASSIFIER_ARTIFACT}"
