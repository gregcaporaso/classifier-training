# classifier-training

Scripts and helpers for training QIIME 2 taxonomic classifiers.

The main entry point is `gtdb-sbatch.sh`, a Slurm batch script that downloads GTDB reference data, trains a naive Bayes classifier, and annotates the resulting artifacts.

## Usage

For step-by-step instructions on configuring and running the script, see the [classifier training procedure](https://github.com/rachis-org/roadmap/blob/main/caporaso/procedures/classifier_training.md).
