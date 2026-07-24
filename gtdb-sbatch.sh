#!/bin/bash
#SBATCH --job-name=train-gtdb
#SBATCH --time=48:00:00
#SBATCH --mem=100G
#SBATCH --partition rhel10

set -e

log_status() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] [train-gtdb] $*"
}

GTDB_VERSION="232.0"
GTDB_DOMAIN="Both"
GTDB_DB_TYPE="SpeciesReps"
GTDB_URL_TYPE="Primary"

# rachis version, matching the value reported as "rachis version" by `qiime info`.
RACHIS_VERSION="$(python -c 'import qiime2; print(qiime2.__version__)')"

GTDB_OUTPUT_PREFIX="gtdb-r${GTDB_VERSION}-${RACHIS_VERSION}"
GTDB_SEQS_ARTIFACT="${GTDB_OUTPUT_PREFIX}-reference-sequences.qza"
GTDB_TAXONOMY_ARTIFACT="${GTDB_OUTPUT_PREFIX}-reference-taxonomy.qza"
# The classifier is always written unsigned; signing is a manual step (see
# below) that renames it to drop the -UNSIGNED suffix. The signed name is
# tracked only so an already-signed classifier still short-circuits a re-run.
GTDB_CLASSIFIER_ARTIFACT="${GTDB_OUTPUT_PREFIX}-classifier-UNSIGNED.qza"
GTDB_CLASSIFIER_SIGNED_ARTIFACT="${GTDB_OUTPUT_PREFIX}-classifier.qza"

# Human-readable description of the reference database for the note annotation,
# derived from the db-type and domain (see `qiime rescript get-gtdb-data --help`).
case "${GTDB_DOMAIN}" in
  Both)     GTDB_DOMAIN_DESCRIPTION="Bacteria and Archaea" ;;
  Bacteria) GTDB_DOMAIN_DESCRIPTION="Bacteria" ;;
  Archaea)  GTDB_DOMAIN_DESCRIPTION="Archaea" ;;
  *)        GTDB_DOMAIN_DESCRIPTION="${GTDB_DOMAIN}" ;;
esac
case "${GTDB_DB_TYPE}" in
  SpeciesReps)
    # 'domain' only applies to the SpeciesReps db-type.
    GTDB_DB_DESCRIPTION="SSU sequences from the set of GTDB representative species for ${GTDB_DOMAIN_DESCRIPTION}" ;;
  All)
    # 'domain' is ignored by GTDB for the 'All' db-type.
    GTDB_DB_DESCRIPTION="all quality-controlled GTDB SSU sequences (not clustered into representative species)" ;;
  *)
    GTDB_DB_DESCRIPTION="GTDB SSU sequences (db-type ${GTDB_DB_TYPE}, domain ${GTDB_DOMAIN_DESCRIPTION})" ;;
esac
GTDB_NOTE_TEXT="GTDB r${GTDB_VERSION} ${GTDB_DB_DESCRIPTION}"

log_status "Pipeline starting for ${GTDB_CLASSIFIER_ARTIFACT}"
log_status "Checking whether classifier artifact already exists"

if [[ -f "${GTDB_CLASSIFIER_ARTIFACT}" || -f "${GTDB_CLASSIFIER_SIGNED_ARTIFACT}" ]]; then
  log_status "Classifier already exists; bailing out without downloading GTDB data or retraining"
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

log_status "Annotating artifacts with reference-database note: ${GTDB_NOTE_TEXT}"
for artifact in \
  "${GTDB_CLASSIFIER_ARTIFACT}" \
  "${GTDB_SEQS_ARTIFACT}" \
  "${GTDB_TAXONOMY_ARTIFACT}"; do
  qiime tools annotation-create \
    --input-path "${artifact}" \
    --annotation-type Note \
    --name reference-database-note \
    --text "${GTDB_NOTE_TEXT}" \
    --output-path "${artifact}"
  log_status "Annotated ${artifact}"
done

# NOTE: the classifier is intentionally left unsigned (hence the -UNSIGNED
# suffix). Signing uses a passphrase-protected GPG key, which cannot be unlocked
# non-interactively in a batch job, so it is done as a manual step afterward
# (see the training procedure docs).
log_status "Pipeline complete: ${GTDB_CLASSIFIER_ARTIFACT} (unsigned; sign manually)"
