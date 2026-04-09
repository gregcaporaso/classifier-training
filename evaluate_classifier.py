#!/usr/bin/env python3

from __future__ import annotations

import html
import json
from datetime import datetime
from pathlib import Path

import click
import numpy as np
import pandas as pd
import qiime2
from qiime2.plugins.feature_classifier.actions import classify_sklearn
from qiime2.plugins.taxa.actions import collapse
from scipy.stats import spearmanr


DEFAULT_MAX_PLOT_POINTS = 10_000
MAX_TAXA_PREVIEW = 50
MAX_FEATURE_DIFF_PREVIEW = 100


def parse_levels(levels_text: str) -> list[int]:
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


def table_artifact_to_taxa_by_sample_df(table_artifact: qiime2.Artifact) -> pd.DataFrame:
    df = table_artifact.view(pd.DataFrame)

    # Collapsed tables should have taxa as rows; transpose when needed.
    index_has_taxonomy_shape = any(";" in str(value) for value in df.index[:10])
    columns_have_taxonomy_shape = any(";" in str(value) for value in df.columns[:10])
    if not index_has_taxonomy_shape and columns_have_taxonomy_shape:
        df = df.transpose()

    df = df.fillna(0)
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)
    return df


def format_taxonomy_for_tooltip(taxon: str) -> str:
    parts = [part.strip() for part in taxon.split(";")]
    parts = [part for part in parts if part]
    if not parts:
        return taxon
    return ";\n".join(parts)


def write_taxa_list_file(path: Path, taxa: list[str]) -> None:
    path.write_text("\n".join(taxa) + "\n", encoding="utf-8")


def make_taxa_preview_html(items: list[str], summary: str, full_filename: str) -> str:
    preview = items[:MAX_TAXA_PREVIEW]
    lines = "\n".join(f"<li>{html.escape(item)}</li>" for item in preview)
    extra = ""
    if len(items) > MAX_TAXA_PREVIEW:
        extra = (
            f"<p>Showing first {MAX_TAXA_PREVIEW} of {len(items)} entries. "
            f"Full list: <a href=\"{html.escape(full_filename)}\">{html.escape(full_filename)}</a></p>"
        )
    else:
        extra = (
            f"<p>Full list: <a href=\"{html.escape(full_filename)}\">{html.escape(full_filename)}</a></p>"
        )

    return (
        f"<details><summary>{html.escape(summary)}</summary>"
        f"{extra}<ul>{lines}</ul></details>"
    )


def build_scatter_plot_df(
    old_shared: pd.DataFrame,
    new_shared: pd.DataFrame,
    shared_taxa: list[str],
    all_samples: list[str],
    max_points: int,
) -> pd.DataFrame:
    old_flat = old_shared.to_numpy(dtype=float).ravel()
    new_flat = new_shared.to_numpy(dtype=float).ravel()
    taxa_labels = np.array([taxon for taxon in shared_taxa for _sample in all_samples], dtype=object)
    sample_labels = np.array(all_samples * len(shared_taxa), dtype=object)

    n_points = len(old_flat)
    if n_points == 0:
        return pd.DataFrame(columns=["x", "y", "taxon_wrapped", "sample", "old_count", "new_count"])

    if n_points > max_points:
        rng = np.random.default_rng(0)
        keep_idx = np.sort(rng.choice(n_points, size=max_points, replace=False))
    else:
        keep_idx = np.arange(n_points)

    old_keep = old_flat[keep_idx]
    new_keep = new_flat[keep_idx]
    taxa_keep = taxa_labels[keep_idx]
    sample_keep = sample_labels[keep_idx]

    return pd.DataFrame(
        {
            "x": np.log1p(old_keep),
            "y": np.log1p(new_keep),
            "taxon_wrapped": [format_taxonomy_for_tooltip(t) for t in taxa_keep],
            "sample": sample_keep,
            "old_count": old_keep,
            "new_count": new_keep,
        }
    )


def write_level_page(
    output_dir: Path,
    level: int,
    old_taxa_count: int,
    new_taxa_count: int,
    shared_taxa_count: int,
    old_total_freq: float,
    new_total_freq: float,
    rho: float,
    p_value: float,
    scatter_plot_records: list[dict[str, object]] | None,
    sampled_points_count: int,
    total_points_count: int,
    scatter_max_value: float,
    shared_taxa: list[str],
    old_only_taxa: list[str],
    new_only_taxa: list[str],
) -> str:
    page_filename = f"level-{level}.html"

    shared_file = f"level-{level}-shared-taxa.txt"
    old_only_file = f"level-{level}-old-only-taxa.txt"
    new_only_file = f"level-{level}-new-only-taxa.txt"

    write_taxa_list_file(output_dir / shared_file, shared_taxa)
    write_taxa_list_file(output_dir / old_only_file, old_only_taxa)
    write_taxa_list_file(output_dir / new_only_file, new_only_taxa)

    if scatter_plot_records is None:
        chart_html = "<p>No shared taxa available for plotting.</p>"
        chart_note = ""
    else:
        chart_spec = {
            "$schema": "https://vega.github.io/schema/vega-lite/v5.json",
            "width": 760,
            "height": 500,
            "layer": [
                {
                    "data": {"values": [{"x": 0.0, "y": 0.0}, {"x": scatter_max_value, "y": scatter_max_value}]},
                    "mark": {"type": "line", "strokeDash": [6, 4], "color": "#555"},
                    "encoding": {
                        "x": {"field": "x", "type": "quantitative"},
                        "y": {"field": "y", "type": "quantitative"},
                    },
                },
                {
                    "data": {"values": scatter_plot_records},
                    "mark": {"type": "circle", "opacity": 0.4, "size": 35},
                    "encoding": {
                        "x": {
                            "field": "x",
                            "type": "quantitative",
                            "title": "Old taxonomy counts per sample-taxon pair (log1p)",
                        },
                        "y": {
                            "field": "y",
                            "type": "quantitative",
                            "title": "New taxonomy counts per sample-taxon pair (log1p)",
                        },
                        "tooltip": [
                            {"field": "taxon_wrapped", "type": "nominal", "title": "Taxon"},
                            {"field": "sample", "type": "nominal", "title": "Sample"},
                            {"field": "old_count", "type": "quantitative", "title": "Old count", "format": ",.0f"},
                            {"field": "new_count", "type": "quantitative", "title": "New count", "format": ",.0f"},
                        ],
                    },
                },
            ],
            "title": f"Level {level} shared taxa frequencies (Spearman rho={rho:.4f}, p-value={p_value:.3e})",
        }

        chart_html = f"""
<div id="level-{level}-scatterplot" class="altair-scatterplot"></div>
<script>
  (function() {{
    const spec = {json.dumps(chart_spec)};
    vegaEmbed('#level-{level}-scatterplot', spec, {{actions: false}}).catch(console.error);
  }})();
</script>
"""
        if sampled_points_count < total_points_count:
            chart_note = (
                f"<p>Plot shows a deterministic sample of {sampled_points_count:,} "
                f"of {total_points_count:,} points to keep page size manageable.</p>"
            )
        else:
            chart_note = f"<p>Plot includes all {total_points_count:,} points.</p>"

    page_html = f"""
<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>Classifier Evaluation - Level {level}</title>
  <script src="https://cdn.jsdelivr.net/npm/vega@5"></script>
  <script src="https://cdn.jsdelivr.net/npm/vega-lite@5"></script>
  <script src="https://cdn.jsdelivr.net/npm/vega-embed@6"></script>
  <style>
    body {{ font-family: -apple-system, BlinkMacSystemFont, Segoe UI, sans-serif; margin: 2rem auto; max-width: 1200px; line-height: 1.5; padding: 0 1rem; }}
    section {{ border: 1px solid #ddd; border-radius: 8px; padding: 1rem 1.25rem; margin: 1rem 0; }}
    h1, h2 {{ margin-top: 0; }}
    .altair-scatterplot {{ border: 1px solid #eee; border-radius: 6px; overflow: hidden; }}
    details {{ margin-top: 0.75rem; }}
    ul {{ max-height: 20rem; overflow-y: auto; background: #fafafa; padding: 0.75rem 1.75rem; border-radius: 6px; border: 1px solid #eee; }}
  </style>
</head>
<body>
  <h1>Taxonomic Level {level}</h1>
  <p><a href="index.html">Back to index</a></p>
  <section>
    <p><strong>Old table taxa:</strong> {old_taxa_count} | <strong>New table taxa:</strong> {new_taxa_count} | <strong>Shared taxa:</strong> {shared_taxa_count}</p>
    <p><strong>Old total frequency:</strong> {old_total_freq:.0f} | <strong>New total frequency:</strong> {new_total_freq:.0f}</p>
    <p><strong>Spearman (shared taxa, sample-wise frequencies):</strong> rho={rho:.4f}, p-value={p_value:.3e}</p>
    {chart_note}
    {chart_html}
    {make_taxa_preview_html(shared_taxa, f'Shared taxa ({len(shared_taxa)})', shared_file)}
    {make_taxa_preview_html(old_only_taxa, f'Only in old taxonomy ({len(old_only_taxa)})', old_only_file)}
    {make_taxa_preview_html(new_only_taxa, f'Only in new taxonomy ({len(new_only_taxa)})', new_only_file)}
  </section>
</body>
</html>
"""

    (output_dir / page_filename).write_text(page_html, encoding="utf-8")
    return page_filename


def format_artifact_header_line(label: str, artifact_fp: Path, artifact: qiime2.Artifact) -> str:
    return (
        f"<p><strong>{html.escape(label)}:</strong> {html.escape(str(artifact_fp))} "
        f"(UUID: {html.escape(str(artifact.uuid))})</p>"
    )


def load_taxonomy_dict(taxonomy_artifact: qiime2.Artifact) -> dict[str, str]:
    """Load FeatureData[Taxonomy] artifact and return dict mapping feature_id to taxonomy."""
    tax_df = taxonomy_artifact.view(pd.DataFrame)
    tax_dict = {}
    for feature_id in tax_df.index:
        if "Taxon" in tax_df.columns:
            tax_str = str(tax_df.loc[feature_id, "Taxon"])
        else:
            tax_str = str(tax_df.loc[feature_id].iloc[0])
        tax_dict[str(feature_id)] = tax_str
    return tax_dict


def render_differing_features_preview_html(differing_features: list[dict[str, str]]) -> str:
    """Render a preview HTML table of features with differing taxonomies."""
    if not differing_features:
        return "<p>No features with differing taxonomies found.</p>"

    preview_features = differing_features[:MAX_FEATURE_DIFF_PREVIEW]
    rows = []
    for feature in preview_features:
        feature_id_escaped = html.escape(feature["feature_id"])
        old_tax_formatted = html.escape(feature["old_taxonomy"]).replace(";", ";<br>")
        new_tax_formatted = html.escape(feature["new_taxonomy"]).replace(";", ";<br>")
        rows.append(
            f"""
  <tr style="border-bottom: 1px solid #eee;">
    <td style="font-family: monospace; font-size: 0.9em; word-break: break-all; padding: 0.5rem; width: 25%;">{feature_id_escaped}</td>
    <td style="padding: 0.5rem; width: 37.5%; border-left: 1px solid #ddd;">{old_tax_formatted}</td>
    <td style="padding: 0.5rem; width: 37.5%; border-left: 1px solid #ddd;">{new_tax_formatted}</td>
  </tr>
"""
        )

    footer = ""
    if len(differing_features) > MAX_FEATURE_DIFF_PREVIEW:
        footer = (
            f"<p>Showing first {MAX_FEATURE_DIFF_PREVIEW} rows. "
            "Use the full TSV for all differing features.</p>"
        )

    table_html = f"""
<div style="overflow-x: auto; max-height: 40rem; overflow-y: auto;">
<table style="width: 100%; border-collapse: collapse; font-size: 0.9em; table-layout: fixed;">
  <thead style="position: sticky; top: 0; background: #fff;">
    <tr style="background: #f5f5f5; border-bottom: 2px solid #ddd;">
      <th style="text-align: left; padding: 0.5rem; font-weight: bold; width: 25%;">Feature ID</th>
      <th style="text-align: left; padding: 0.5rem; font-weight: bold; width: 37.5%; border-left: 1px solid #ddd;">Old Taxonomy</th>
      <th style="text-align: left; padding: 0.5rem; font-weight: bold; width: 37.5%; border-left: 1px solid #ddd;">New Taxonomy</th>
    </tr>
  </thead>
  <tbody>
    {''.join(rows)}
  </tbody>
</table>
</div>
{footer}
"""
    return table_html


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--classifier",
    "classifier_fp",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Path to trained classifier artifact (.qza).",
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
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Path to FeatureData[Sequence] artifact (.qza).",
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
    default=DEFAULT_MAX_PLOT_POINTS,
    show_default=True,
    type=click.IntRange(min=1),
    help="Maximum number of points to include per level scatter plot.",
)
def main(
    classifier_fp: Path,
    old_taxonomy_fp: Path,
    sequences_fp: Path,
    table_fp: Path,
    output_dir: Path,
    levels: str,
    confidence: float,
    max_plot_points: int,
) -> None:
    """Evaluate a newly trained classifier against an older taxonomy assignment."""
    comparison_levels = parse_levels(levels)
    output_dir.mkdir(parents=True, exist_ok=True)

    click.echo("Loading artifacts...")
    classifier = qiime2.Artifact.load(str(classifier_fp))
    old_taxonomy = qiime2.Artifact.load(str(old_taxonomy_fp))
    sequences = qiime2.Artifact.load(str(sequences_fp))
    table = qiime2.Artifact.load(str(table_fp))

    click.echo("Classifying sequences with q2-feature-classifier classify-sklearn...")
    new_taxonomy, = classify_sklearn(
        reads=sequences,
        classifier=classifier,
        confidence=confidence,
    )
    new_taxonomy_fp = output_dir / "new-taxonomy.qza"
    new_taxonomy.save(str(new_taxonomy_fp))

    click.echo("Comparing feature-level taxonomies...")
    old_tax_dict = load_taxonomy_dict(old_taxonomy)
    new_tax_dict = load_taxonomy_dict(new_taxonomy)

    all_features = set(old_tax_dict.keys()) | set(new_tax_dict.keys())
    differing_features = []
    for feature_id in sorted(all_features):
        old_tax = old_tax_dict.get(feature_id, "")
        new_tax = new_tax_dict.get(feature_id, "")
        if old_tax != new_tax:
            differing_features.append({
                "feature_id": feature_id,
                "old_taxonomy": old_tax,
                "new_taxonomy": new_tax,
            })

    level_index_rows: list[str] = []

    for level in comparison_levels:
        click.echo(f"Collapsing and comparing level {level}...")

        old_collapsed, = collapse(table=table, taxonomy=old_taxonomy, level=level)
        new_collapsed, = collapse(table=table, taxonomy=new_taxonomy, level=level)

        old_collapsed_fp = output_dir / f"old-collapsed-l{level}.qza"
        new_collapsed_fp = output_dir / f"new-collapsed-l{level}.qza"
        old_collapsed.save(str(old_collapsed_fp))
        new_collapsed.save(str(new_collapsed_fp))

        old_df = table_artifact_to_taxa_by_sample_df(old_collapsed)
        new_df = table_artifact_to_taxa_by_sample_df(new_collapsed)

        old_taxa = set(old_df.index)
        new_taxa = set(new_df.index)
        shared_taxa = sorted(old_taxa & new_taxa)
        old_only_taxa = sorted(old_taxa - new_taxa)
        new_only_taxa = sorted(new_taxa - old_taxa)

        all_samples = sorted(set(old_df.columns) | set(new_df.columns))
        old_shared = old_df.reindex(index=shared_taxa, columns=all_samples, fill_value=0)
        new_shared = new_df.reindex(index=shared_taxa, columns=all_samples, fill_value=0)

        old_flat = old_shared.to_numpy(dtype=float).ravel()
        new_flat = new_shared.to_numpy(dtype=float).ravel()
        if len(old_flat) > 1:
            rho, p_value = spearmanr(old_flat, new_flat)
        else:
            rho, p_value = float("nan"), float("nan")

        if len(old_flat) > 0 and len(shared_taxa) > 0:
            plot_df = build_scatter_plot_df(
                old_shared=old_shared,
                new_shared=new_shared,
                shared_taxa=shared_taxa,
                all_samples=all_samples,
                max_points=max_plot_points,
            )
            scatter_plot_records = plot_df.to_dict(orient="records")
            sampled_points_count = len(plot_df)
            total_points_count = len(old_flat)
            if len(plot_df) > 0:
                scatter_max_value = float(max(plot_df["x"].max(), plot_df["y"].max()))
            else:
                scatter_max_value = 0.0
        else:
            scatter_plot_records = None
            sampled_points_count = 0
            total_points_count = len(old_flat)
            scatter_max_value = 0.0

        level_page_filename = write_level_page(
            output_dir=output_dir,
            level=level,
            old_taxa_count=len(old_taxa),
            new_taxa_count=len(new_taxa),
            shared_taxa_count=len(shared_taxa),
            old_total_freq=float(old_df.to_numpy(dtype=float).sum()),
            new_total_freq=float(new_df.to_numpy(dtype=float).sum()),
            rho=rho,
            p_value=p_value,
            scatter_plot_records=scatter_plot_records,
            sampled_points_count=sampled_points_count,
            total_points_count=total_points_count,
            scatter_max_value=scatter_max_value,
            shared_taxa=shared_taxa,
            old_only_taxa=old_only_taxa,
            new_only_taxa=new_only_taxa,
        )

        level_index_rows.append(
            f"<tr><td><a href=\"{level_page_filename}\">Level {level}</a></td>"
            f"<td>{len(old_taxa)}</td><td>{len(new_taxa)}</td><td>{len(shared_taxa)}</td>"
            f"<td>{rho:.4f}</td><td>{p_value:.3e}</td></tr>"
        )

    feature_diff_df = pd.DataFrame(differing_features)
    feature_diff_tsv = output_dir / "feature-taxonomy-differences.tsv"
    if len(feature_diff_df) == 0:
        feature_diff_df = pd.DataFrame(columns=["feature_id", "old_taxonomy", "new_taxonomy"])
    feature_diff_df.to_csv(feature_diff_tsv, sep="\t", index=False)

    feature_comparison_preview = render_differing_features_preview_html(differing_features)

    index_html = f"""
<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
    <title>Classifier Evaluation Report Index</title>
  <style>
    body {{ font-family: -apple-system, BlinkMacSystemFont, Segoe UI, sans-serif; margin: 2rem auto; max-width: 1200px; line-height: 1.5; padding: 0 1rem; }}
    section {{ border: 1px solid #ddd; border-radius: 8px; padding: 1rem 1.25rem; margin: 1rem 0; }}
    h1, h2 {{ margin-top: 0; }}
        table {{ width: 100%; border-collapse: collapse; }}
        th, td {{ border-bottom: 1px solid #eee; padding: 0.5rem; text-align: left; }}
        thead tr {{ background: #f5f5f5; }}
    details {{ margin-top: 0.75rem; }}
    ul {{ max-height: 20rem; overflow-y: auto; background: #fafafa; padding: 0.75rem 1.75rem; border-radius: 6px; border: 1px solid #eee; }}
  </style>
</head>
<body>
    <h1>Classifier Evaluation Report</h1>
  <p><strong>Generated:</strong> {html.escape(datetime.now().isoformat(timespec='seconds'))}</p>
    {format_artifact_header_line('Classifier', classifier_fp, classifier)}
    {format_artifact_header_line('Old taxonomy', old_taxonomy_fp, old_taxonomy)}
    {format_artifact_header_line('Sequences', sequences_fp, sequences)}
    {format_artifact_header_line('Table', table_fp, table)}
  <p><strong>Levels:</strong> {', '.join(map(str, comparison_levels))}</p>
    {format_artifact_header_line('New taxonomy output', new_taxonomy_fp, new_taxonomy)}

    <section>
        <h2>Level Comparisons</h2>
        <table>
            <thead>
                <tr>
                    <th>Level Page</th>
                    <th>Old taxa</th>
                    <th>New taxa</th>
                    <th>Shared taxa</th>
                    <th>Spearman rho</th>
                    <th>p-value</th>
                </tr>
            </thead>
            <tbody>
                {''.join(level_index_rows)}
            </tbody>
        </table>
    </section>

    <section>
        <h2>Feature-level Taxonomy Comparison</h2>
        <p><strong>Total features:</strong> {len(all_features)} | <strong>Features with differing taxonomy:</strong> {len(differing_features)}</p>
        <p><strong>Full table:</strong> <a href="feature-taxonomy-differences.tsv">feature-taxonomy-differences.tsv</a></p>
        {feature_comparison_preview}
    </section>
</body>
</html>
"""

    index_fp = output_dir / "index.html"
    index_fp.write_text(index_html, encoding="utf-8")
    click.echo(f"Report bundle index written: {index_fp}")


if __name__ == "__main__":
    main()
