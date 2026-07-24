#!/usr/bin/env python3

from __future__ import annotations

import subprocess
import tempfile
import uuid
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import click
import pandas as pd
import qiime2
from qiime2.plugins.feature_classifier.actions import classify_sklearn
from qiime2.plugins.metadata.visualizers import tabulate  # type: ignore[import]
from qiime2.plugins.taxa.actions import collapse
from qiime2.plugins.vizard.actions import scatterplot_2d  # type: ignore[import]

from _actions import build_scatter_metadata, compare_collapsed_tables, compare_taxonomy


def _parse_levels(levels_text: str) -> list[int]:
    levels: list[int] = []
    for token in levels_text.split(","):
        token = token.strip()
        if not token:
            continue
        try:
            level = int(token)
        except ValueError as e:
            raise click.BadParameter(
                f"Invalid taxonomic level '{token}'. Levels must be integers."
            ) from e
        if level < 1:
            raise click.BadParameter("Taxonomic levels must be >= 1.")
        levels.append(level)

    if not levels:
        raise click.BadParameter("At least one level is required.")

    return sorted(set(levels))


def _write_visualization(viz_dir: Path, output_path: Path) -> None:
    """Package a directory of HTML/asset files as a minimal valid QIIME 2 .qzv.

    This replicates the QIIME 2 archive v5 format so the resulting file can be
    loaded by ``qiime tools make-report`` and QIIME 2 View.  Once
    ``compare_collapsed_tables`` is registered as a plugin visualizer this
    helper is no longer needed — the framework will handle packaging.
    """
    viz_uuid = str(uuid.uuid4())
    exec_uuid = str(uuid.uuid4())
    now = datetime.now(timezone.utc).isoformat(timespec="seconds")

    with zipfile.ZipFile(output_path, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr(
            f"{viz_uuid}/VERSION",
            "QIIME 2\narchive: 5\nframework: 2024.10.0\n",
        )
        zf.writestr(
            f"{viz_uuid}/metadata.yaml",
            f"uuid: {viz_uuid}\ntype: Visualization\nformat: null\n",
        )
        zf.writestr(f"{viz_uuid}/provenance/citations.bib", "")
        zf.writestr(
            f"{viz_uuid}/provenance/metadata.yaml",
            (
                f"execution:\n"
                f"  uuid: {exec_uuid}\n"
                f"  runtime:\n"
                f"    start: '{now}'\n"
                f"    end: '{now}'\n"
                f"    duration: PT0S\n"
                f"artifact:\n"
                f"  uuid: {viz_uuid}\n"
                f"  type: Visualization\n"
                f"  format: null\n"
            ),
        )
        zf.writestr(
            f"{viz_uuid}/provenance/action/action.yaml",
            "action:\n  type: visualizer\noutput-name: visualization\ncitations: []\n",
        )
        for f in sorted(viz_dir.rglob("*")):
            if f.is_file():
                zf.write(f, f"{viz_uuid}/data/{f.relative_to(viz_dir)}")


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--classifier",
    "classifier_fp",
    default=None,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Path to trained classifier artifact (.qza). "
         "Mutually exclusive with --new-taxonomy.",
)
@click.option(
    "--new-taxonomy",
    "new_taxonomy_fp",
    default=None,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Path to a pre-computed FeatureData[Taxonomy] artifact (.qza) to use "
         "instead of running classify-sklearn. "
         "Mutually exclusive with --classifier.",
)
@click.option(
    "--old-taxonomy",
    "old_taxonomy_fp",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Path to old taxonomy artifact FeatureData[Taxonomy] (.qza).",
)
@click.option(
    "--sequences",
    "sequences_fp",
    default=None,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Path to FeatureData[Sequence] artifact (.qza). "
         "Required when --classifier is used.",
)
@click.option(
    "--table",
    "table_fp",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Path to FeatureTable[Frequency] artifact (.qza).",
)
@click.option(
    "--output-dir",
    "output_dir",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
    help="Directory where outputs and report are written.",
)
@click.option(
    "--levels",
    default="2,3,4,5,6,7",
    show_default=True,
    help="Comma-separated taxonomic levels to collapse and compare.",
)
@click.option(
    "--confidence",
    default=0.7,
    show_default=True,
    type=float,
    help="Confidence value passed to classify-sklearn.",
)
@click.option(
    "--max-plot-points",
    default=10_000,
    show_default=True,
    type=click.IntRange(min=1),
    help="Maximum scatter-plot points per level comparison.",
)
def main(
    classifier_fp: Path | None,
    new_taxonomy_fp: Path | None,
    old_taxonomy_fp: Path,
    sequences_fp: Path | None,
    table_fp: Path,
    output_dir: Path,
    levels: str,
    confidence: float,
    max_plot_points: int,
) -> None:
    """Evaluate a newly trained classifier against an older taxonomy assignment."""
    if classifier_fp is None and new_taxonomy_fp is None:
        raise click.UsageError("Provide either --classifier or --new-taxonomy.")
    if classifier_fp is not None and new_taxonomy_fp is not None:
        raise click.UsageError("--classifier and --new-taxonomy are mutually exclusive.")
    if classifier_fp is not None and sequences_fp is None:
        raise click.UsageError("--sequences is required when --classifier is used.")

    comparison_levels = _parse_levels(levels)
    output_dir.mkdir(parents=True, exist_ok=True)

    click.echo("Loading artifacts...")
    old_taxonomy = qiime2.Artifact.load(str(old_taxonomy_fp))
    table = qiime2.Artifact.load(str(table_fp))

    if new_taxonomy_fp is not None:
        click.echo("Loading pre-computed new taxonomy...")
        new_taxonomy = qiime2.Artifact.load(str(new_taxonomy_fp))
    else:
        sequences = qiime2.Artifact.load(str(sequences_fp))
        classifier = qiime2.Artifact.load(str(classifier_fp))
        click.echo("Classifying sequences with classify-sklearn...")
        (new_taxonomy,) = classify_sklearn(
            reads=sequences,
            classifier=classifier,
            confidence=confidence,
        )
        new_taxonomy.save(str(output_dir / "new-taxonomy.qza"))

    click.echo("Comparing feature-level taxonomies...")
    old_tax_df = old_taxonomy.view(pd.DataFrame)
    new_tax_df = new_taxonomy.view(pd.DataFrame)
    comparison_df = compare_taxonomy(old_tax_df, new_tax_df)
    click.echo(f"  {len(comparison_df)} features with differing taxonomy assignments.")
    comparison_metadata = qiime2.Metadata(comparison_df)
    (taxonomy_comparison_viz,) = tabulate(input=comparison_metadata)
    taxonomy_comparison_viz.save(str(output_dir / "taxonomy-comparison.qzv"))

    qzv_paths: list[Path] = [output_dir / "taxonomy-comparison.qzv"]

    for level in comparison_levels:
        click.echo(f"Collapsing and comparing level {level}...")

        (old_collapsed,) = collapse(table=table, taxonomy=old_taxonomy, level=level)
        (new_collapsed,) = collapse(table=table, taxonomy=new_taxonomy, level=level)
        old_collapsed.save(str(output_dir / f"old-collapsed-l{level}.qza"))
        new_collapsed.save(str(output_dir / f"new-collapsed-l{level}.qza"))

        old_df = old_collapsed.view(pd.DataFrame)
        new_df = new_collapsed.view(pd.DataFrame)

        with tempfile.TemporaryDirectory() as tmp:
            compare_collapsed_tables(
                output_dir=tmp,
                old_table=old_df,
                new_table=new_df,
                level=level,
            )
            stats_qzv_path = output_dir / f"level-{level}-comparison.qzv"
            _write_visualization(Path(tmp), stats_qzv_path)

        scatter_df = build_scatter_metadata(old_df, new_df, max_plot_points)
        scatter_metadata = qiime2.Metadata(scatter_df)
        (scatter_viz,) = scatterplot_2d(
            metadata=scatter_metadata,
            x_measure="log1p-old-count",
            y_measure="log1p-new-count",
            title=f"Level {level} shared-taxa frequencies (log1p)",
        )
        scatter_qzv_path = output_dir / f"level-{level}-scatter.qzv"
        scatter_viz.save(str(scatter_qzv_path))

        qzv_paths.extend([stats_qzv_path, scatter_qzv_path])

    click.echo("Building combined report with qiime tools make-report...")
    report_dir = output_dir / "report"
    cmd = ["qiime", "tools", "make-report", "--report-path", str(report_dir)]
    for qzv in qzv_paths:
        cmd += [str(qzv)]
    subprocess.run(cmd, check=True)
    click.echo(f"Report written to: {report_dir}")


if __name__ == "__main__":
    main()
