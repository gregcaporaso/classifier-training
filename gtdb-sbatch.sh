#!/bin/bash
#SBATCH --job-name=train-gtdb
#SBATCH --time=48:00:00
#SBATCH --mem=100G

set -e

GTDB_VERSION=226
GTDB_BASE_URL="https://data.gtdb.ecogenomic.org/releases/release${GTDB_VERSION}/${GTDB_VERSION}.0"

GTDB_FILES=(
  "genomic_files_reps/ar53_ssu_reps_r${GTDB_VERSION}.fna.gz"
  "genomic_files_reps/bac120_ssu_reps_r${GTDB_VERSION}.fna.gz"
  "bac120_taxonomy_r${GTDB_VERSION}.tsv.gz"
  "ar53_taxonomy_r${GTDB_VERSION}.tsv.gz"
)

for rel_path in "${GTDB_FILES[@]}"; do
  gz_file="$(basename "${rel_path}")"
  input_file="${gz_file%.gz}"

  if [[ ! -f "${input_file}" ]]; then
    if [[ ! -f "${gz_file}" ]]; then
      wget "${GTDB_BASE_URL}/${rel_path}" -O "${gz_file}"
    fi
    gunzip -c "${gz_file}" > "${input_file}"
  fi
done


qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path ar53_ssu_reps_r${GTDB_VERSION}.fna \
  --output-path ar53_ssu_reps.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path ar53_taxonomy_r${GTDB_VERSION}.tsv \
  --output-path ar53_taxonomy.qza

qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path  bac120_ssu_reps_r${GTDB_VERSION}.fna \
  --output-path bac120_ssu_reps.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path bac120_taxonomy_r${GTDB_VERSION}.tsv \
  --output-path bac120_taxonomy.qza

qiime feature-table merge-seqs \
  --i-data ar53_ssu_reps.qza bac120_ssu_reps.qza \
  --o-merged-data gtdb_seqs.qza

qiime feature-table merge-taxa \
   --i-data ar53_taxonomy.qza bac120_taxonomy.qza \
   --o-merged-data gtdb_taxonomy.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads gtdb_seqs.qza \
  --i-reference-taxonomy gtdb_taxonomy.qza \
  --o-classifier gtdb_classifier_r${GTDB_VERSION}.qza
