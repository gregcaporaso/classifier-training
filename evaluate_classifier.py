#!/usr/bin/env python3

from __future__ import annotations

import html
import json
from datetime import datetime
from pathlib import Path

import altair as alt
import click
import numpy as np
import pandas as pd
import qiime2
from qiime2.plugins.feature_classifier.actions import classify_sklearn
from qiime2.plugins.taxa.actions import collapse
from scipy.stats import spearmanr


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


def render_interactive_scatterplot_html(
    old_values: np.ndarray,
    new_values: np.ndarray,
    taxa_labels: list[str],
    sample_labels: list[str],
    level: int,
    rho: float,
    p_value: float,
) -> str:
    x_plot = np.log1p(old_values)
    y_plot = np.log1p(new_values)
    plot_df = pd.DataFrame(
        {
            "x": x_plot,
            "y": y_plot,
            "taxon": taxa_labels,
            "taxon_wrapped": [format_taxonomy_for_tooltip(t) for t in taxa_labels],
            "sample": sample_labels,
            "old_count": old_values,
            "new_count": new_values,
        }
    )

    max_value = float(max(np.max(x_plot, initial=0.0), np.max(y_plot, initial=0.0)))
    diag_df = pd.DataFrame({"x": [0.0, max_value], "y": [0.0, max_value]})

    points = (
        alt.Chart(plot_df)
        .mark_circle(size=35, opacity=0.4)
        .encode(
            x=alt.X("x:Q", title="Old taxonomy counts per sample-taxon pair (log1p)"),
            y=alt.Y("y:Q", title="New taxonomy counts per sample-taxon pair (log1p)"),
            tooltip=[
                alt.Tooltip("taxon_wrapped:N", title="Taxon"),
                alt.Tooltip("sample:N", title="Sample"),
                alt.Tooltip("old_count:Q", title="Old count", format=",.0f"),
                alt.Tooltip("new_count:Q", title="New count", format=",.0f"),
            ],
        )
    )

    diagonal = (
        alt.Chart(diag_df)
        .mark_line(strokeDash=[6, 4], color="#555")
        .encode(x="x:Q", y="y:Q")
    )

    chart = (
        (diagonal + points)
        .properties(
            width=760,
            height=500,
            title=(
                f"Level {level} shared taxa frequencies "
                f"(Spearman rho={rho:.4f}, p-value={p_value:.3e})"
            ),
        )
        .interactive()
    )

    # Large tables are common here; disable the default row cap for Vega-Lite export.
    alt.data_transformers.disable_max_rows()
    chart_spec = json.dumps(chart.to_dict())
    chart_id = f"level-{level}-scatterplot"
    return f"""
<div id="{chart_id}" class="altair-scatterplot"></div>
<script>
    (function() {{
        const spec = {chart_spec};
        vegaEmbed('#{chart_id}', spec, {{actions: false}}).catch(console.error);
    }})();
</script>
"""


def make_taxa_list_html(items: list[str], summary: str) -> str:
    lines = "\n".join(f"<li>{html.escape(item)}</li>" for item in items)
    return (
        f"<details><summary>{html.escape(summary)}</summary>"
        f"<ul>{lines}</ul>"
        f"</details>"
    )


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
        tax_str = str(tax_df.loc[feature_id, "Taxon"]) if "Taxon" in tax_df.columns else str(tax_df.iloc[feature_id, 0])
        tax_dict[str(feature_id)] = tax_str
    return tax_dict


def render_differing_features_table_html(differing_features: list[dict[str, str]]) -> str:
    """Render an HTML table of features with differing taxonomies."""
    if not differing_features:
        return "<p>No features with differing taxonomies found.</p>"

    rows = []
    for feature in differing_features:
        feature_id_escaped = html.escape(feature["feature_id"])
        old_tax_formatted = feature["old_taxonomy"].replace(";", ";<br>")
        new_tax_formatted = feature["new_taxonomy"].replace(";", ";<br>")
        rows.append(
            f"""
  <tr style="border-bottom: 1px solid #eee;">
    <td style="font-family: monospace; font-size: 0.9em; word-break: break-all; padding: 0.5rem; width: 25%;">{feature_id_escaped}</td>
    <td style="padding: 0.5rem; width: 37.5%; border-left: 1px solid #ddd;">{old_tax_formatted}</td>
    <td style="padding: 0.5rem; width: 37.5%; border-left: 1px solid #ddd;">{new_tax_formatted}</td>
  </tr>
"""
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
def main(
    classifier_fp: Path,
    old_taxonomy_fp: Path,
    sequences_fp: Path,
    table_fp: Path,
    output_dir: Path,
    levels: str,
    confidence: float,
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

    report_sections: list[str] = []

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
        taxa_labels = [taxon for taxon in shared_taxa for _sample in all_samples]
        sample_labels = all_samples * len(shared_taxa)
        if len(old_flat) > 1:
            rho, p_value = spearmanr(old_flat, new_flat)
        else:
            rho, p_value = float("nan"), float("nan")

        if len(old_flat) > 0:
            plot_html = render_interactive_scatterplot_html(
                old_values=old_flat,
                new_values=new_flat,
                taxa_labels=taxa_labels,
                sample_labels=sample_labels,
                level=level,
                rho=rho,
                p_value=p_value,
            )
        else:
            plot_html = "<p>No shared taxa available for plotting.</p>"

        level_section = f"""
<section>
  <h2>Taxonomic level {level}</h2>
  <p><strong>Old table taxa:</strong> {len(old_taxa)} | <strong>New table taxa:</strong> {len(new_taxa)} | <strong>Shared taxa:</strong> {len(shared_taxa)}</p>
  <p><strong>Old total frequency:</strong> {old_df.to_numpy(dtype=float).sum():.0f} | <strong>New total frequency:</strong> {new_df.to_numpy(dtype=float).sum():.0f}</p>
  <p><strong>Spearman (shared taxa, sample-wise frequencies):</strong> rho={rho:.4f}, p-value={p_value:.3e}</p>
  {plot_html}
  {make_taxa_list_html(shared_taxa, f'Shared taxa ({len(shared_taxa)})')}
  {make_taxa_list_html(old_only_taxa, f'Only in old taxonomy ({len(old_only_taxa)})')}
  {make_taxa_list_html(new_only_taxa, f'Only in new taxonomy ({len(new_only_taxa)})')}
</section>
"""
        report_sections.append(level_section)

    feature_comparison_html = render_differing_features_table_html(differing_features)
    feature_comparison_section = f"""
<section>
  <h2>Feature-level Taxonomy Comparison</h2>
  <p><strong>Total features:</strong> {len(all_features)} | <strong>Features with differing taxonomy:</strong> {len(differing_features)}</p>
  {feature_comparison_html}
</section>
"""
    report_sections.append(feature_comparison_section)

    report_html = f"""
<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>Classifier Evaluation Report</title>
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
  <h1>Classifier Evaluation Report</h1>
  <p><strong>Generated:</strong> {html.escape(datetime.now().isoformat(timespec='seconds'))}</p>
    {format_artifact_header_line('Classifier', classifier_fp, classifier)}
    {format_artifact_header_line('Old taxonomy', old_taxonomy_fp, old_taxonomy)}
    {format_artifact_header_line('Sequences', sequences_fp, sequences)}
    {format_artifact_header_line('Table', table_fp, table)}
  <p><strong>Levels:</strong> {', '.join(map(str, comparison_levels))}</p>
    {format_artifact_header_line('New taxonomy output', new_taxonomy_fp, new_taxonomy)}
  {''.join(report_sections)}
</body>
</html>
"""

    report_fp = output_dir / "classifier-evaluation-report.html"
    report_fp.write_text(report_html, encoding="utf-8")
    click.echo(f"Report written: {report_fp}")


if __name__ == "__main__":
    main()
