"""
Functions ready to be registered as QIIME 2 plugin actions.

    compare_taxonomy         → method     → ImmutableMetadata
    compare_collapsed_tables → visualizer → writes index.html
    build_scatter_metadata   → helper for calling scatterplot_2d per level

Registration notes are in each function's docstring.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


# ---------------------------------------------------------------------------
# compare_taxonomy  (method)
# ---------------------------------------------------------------------------

def compare_taxonomy(
    old_taxonomy: pd.DataFrame,
    new_taxonomy: pd.DataFrame,
) -> pd.DataFrame:
    """Compare two FeatureData[Taxonomy] inputs at the feature level.

    Returns a metadata-compatible DataFrame containing **only** features whose
    taxonomy assignments differ.  Columns:

        old-taxon              full old taxonomy string
        new-taxon              full new taxonomy string
        first-differing-level  1-based level index of the first differing token

    Because only differing features are present, passing this artifact to
    ``qiime metadata tabulate`` will — via inner-join merging — automatically
    restrict the table to those features.

    Registration sketch::

        plugin.methods.register_function(
            function=compare_taxonomy,
            inputs={
                'old_taxonomy': FeatureData[Taxonomy],
                'new_taxonomy': FeatureData[Taxonomy],
            },
            parameters={},
            outputs=[('taxonomy_comparison', ImmutableMetadata)],
            input_descriptions={...},
            output_descriptions={...},
            name='Compare taxonomy assignments',
            description='...',
        )
    """

    def _taxon_series(df: pd.DataFrame) -> pd.Series:
        col = 'Taxon' if 'Taxon' in df.columns else df.columns[0]
        return df[col].astype(str)

    old_s = _taxon_series(old_taxonomy)
    new_s = _taxon_series(new_taxonomy)

    common_ids = old_s.index.intersection(new_s.index)
    records = []

    for fid in common_ids:
        old_tax = old_s[fid]
        new_tax = new_s[fid]
        if old_tax == new_tax:
            continue

        old_parts = [p.strip() for p in old_tax.split(';')]
        new_parts = [p.strip() for p in new_tax.split(';')]

        first_diff = None
        for i, (a, b) in enumerate(zip(old_parts, new_parts)):
            if a != b:
                first_diff = i + 1
                break
        if first_diff is None:
            first_diff = min(len(old_parts), len(new_parts)) + 1

        records.append({
            'id': fid,
            'old-taxon': old_tax,
            'new-taxon': new_tax,
            'first-differing-level': first_diff,
        })

    if not records:
        result = pd.DataFrame(columns=['old-taxon', 'new-taxon', 'first-differing-level'])
        result.index.name = 'feature-id'
        return result

    result = pd.DataFrame(records).set_index('id')
    result.index.name = 'feature-id'
    return result


# ---------------------------------------------------------------------------
# compare_collapsed_tables  (visualizer)
# ---------------------------------------------------------------------------

def compare_collapsed_tables(
    output_dir: str,
    old_table: pd.DataFrame,
    new_table: pd.DataFrame,
    level: int,
) -> None:
    """Visualize agreement between two collapsed FeatureTable[Frequency] artifacts.

    Both tables should have already been collapsed to ``level`` via
    ``qiime taxa collapse``.  Reports Spearman correlation, taxa counts, and
    per-category taxa lists.  Use ``build_scatter_metadata`` + ``scatterplot_2d``
    separately to generate the scatter plot visualization.

    Registration sketch::

        plugin.visualizers.register_function(
            function=compare_collapsed_tables,
            inputs={
                'old_table': FeatureTable[Frequency],
                'new_table': FeatureTable[Frequency],
            },
            parameters={
                'level': Int % Range(1, None),
            },
            input_descriptions={
                'old_table': 'Collapsed table from the reference taxonomy.',
                'new_table': 'Collapsed table from the new classifier taxonomy.',
            },
            parameter_descriptions={
                'level': 'Taxonomic level at which both tables were collapsed.',
            },
            name='Compare collapsed feature tables',
            description='...',
        )
    """
    old_df = _orient_taxa_as_rows(old_table)
    new_df = _orient_taxa_as_rows(new_table)

    old_taxa = set(old_df.index)
    new_taxa = set(new_df.index)
    shared_taxa = sorted(old_taxa & new_taxa)
    old_only = sorted(old_taxa - new_taxa)
    new_only = sorted(new_taxa - old_taxa)

    all_samples = sorted(set(old_df.columns) | set(new_df.columns))
    old_shared = old_df.reindex(index=shared_taxa, columns=all_samples, fill_value=0)
    new_shared = new_df.reindex(index=shared_taxa, columns=all_samples, fill_value=0)

    old_flat = old_shared.to_numpy(dtype=float).ravel()
    new_flat = new_shared.to_numpy(dtype=float).ravel()

    if len(old_flat) > 1:
        rho, p_value = spearmanr(old_flat, new_flat)
    else:
        rho, p_value = float('nan'), float('nan')

    out = Path(output_dir)

    for suffix, items in (
        ('shared', shared_taxa),
        ('old-only', old_only),
        ('new-only', new_only),
    ):
        (out / f'level-{level}-{suffix}-taxa.txt').write_text(
            '\n'.join(items) + '\n', encoding='utf-8'
        )

    index_html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8"/>
  <meta name="viewport" content="width=device-width, initial-scale=1"/>
  <title>Collapsed-table comparison – Level {level}</title>
  <style>
    body {{font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
          margin: 2rem auto; max-width: 1100px; line-height: 1.5; padding: 0 1rem;}}
    section {{border: 1px solid #ddd; border-radius: 8px;
              padding: 1rem 1.25rem; margin: 1rem 0;}}
    h1 {{margin-top: 0;}}
    details {{margin-top: 0.75rem;}}
    ul {{max-height: 18rem; overflow-y: auto; background: #fafafa;
         padding: 0.75rem 1.75rem; border-radius: 6px; border: 1px solid #eee;}}
  </style>
</head>
<body>
  <h1>Level {level} – collapsed-table comparison</h1>
  <section>
    <p>
      <strong>Old taxa:</strong> {len(old_taxa)} &nbsp;|&nbsp;
      <strong>New taxa:</strong> {len(new_taxa)} &nbsp;|&nbsp;
      <strong>Shared:</strong> {len(shared_taxa)}
    </p>
    <p>
      <strong>Old total frequency:</strong> {old_flat.sum():,.0f} &nbsp;|&nbsp;
      <strong>New total frequency:</strong> {new_flat.sum():,.0f}
    </p>
    <p>
      <strong>Spearman rho (shared taxa, sample-wise):</strong>
      {rho:.4f} &nbsp; <strong>p-value:</strong> {p_value:.3e}
    </p>
    {_taxa_details(shared_taxa, f'Shared taxa ({len(shared_taxa)})',
                   f'level-{level}-shared-taxa.txt')}
    {_taxa_details(old_only, f'Old-only taxa ({len(old_only)})',
                   f'level-{level}-old-only-taxa.txt')}
    {_taxa_details(new_only, f'New-only taxa ({len(new_only)})',
                   f'level-{level}-new-only-taxa.txt')}
  </section>
</body>
</html>
"""
    (out / 'index.html').write_text(index_html, encoding='utf-8')


# ---------------------------------------------------------------------------
# build_scatter_metadata  (helper — not a plugin action)
# ---------------------------------------------------------------------------

_MAX_SCATTER_POINTS = 10_000


def build_scatter_metadata(
    old_table: pd.DataFrame,
    new_table: pd.DataFrame,
    max_points: int = _MAX_SCATTER_POINTS,
) -> pd.DataFrame:
    """Build a metadata-compatible DataFrame for use with scatterplot_2d.

    Each row is a (sample, taxon) pair for shared taxa only.  Columns:

        log1p-old-count   log1p-transformed old taxonomy frequency
        log1p-new-count   log1p-transformed new taxonomy frequency

    The index (``id``) is ``"{sample}|{taxon}"``.  If the number of pairs
    exceeds ``max_points`` a deterministic random subsample is taken.
    """
    old_df = _orient_taxa_as_rows(old_table)
    new_df = _orient_taxa_as_rows(new_table)

    shared_taxa = sorted(set(old_df.index) & set(new_df.index))
    if not shared_taxa:
        result = pd.DataFrame(columns=['log1p-old-count', 'log1p-new-count'])
        result.index.name = 'id'
        return result

    all_samples = sorted(set(old_df.columns) | set(new_df.columns))
    old_shared = old_df.reindex(index=shared_taxa, columns=all_samples, fill_value=0)
    new_shared = new_df.reindex(index=shared_taxa, columns=all_samples, fill_value=0)

    old_flat = old_shared.to_numpy(dtype=float).ravel()
    new_flat = new_shared.to_numpy(dtype=float).ravel()
    ids = np.array(
        [f'{sample}|{taxon}' for taxon in shared_taxa for sample in all_samples],
        dtype=object,
    )

    n = len(old_flat)
    if n > max_points:
        rng = np.random.default_rng(0)
        idx = np.sort(rng.choice(n, size=max_points, replace=False))
        old_flat, new_flat, ids = old_flat[idx], new_flat[idx], ids[idx]

    result = pd.DataFrame(
        {
            'log1p-old-count': np.log1p(old_flat),
            'log1p-new-count': np.log1p(new_flat),
        },
        index=pd.Index(ids, name='id'),
    )
    return result


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _orient_taxa_as_rows(df: pd.DataFrame) -> pd.DataFrame:
    """Ensure taxa strings are the row index (transpose if needed)."""
    index_has_taxa = any(';' in str(v) for v in df.index[:10])
    cols_have_taxa = any(';' in str(v) for v in df.columns[:10])
    if not index_has_taxa and cols_have_taxa:
        df = df.T
    return df.fillna(0).astype(float)


def _taxa_details(items: list[str], summary: str, filename: str, preview: int = 50) -> str:
    shown = items[:preview]
    rows = '\n'.join(f'<li>{item}</li>' for item in shown)
    more = (
        f'<p>Showing {preview} of {len(items)}. '
        f'Full list: <a href="{filename}">{filename}</a></p>'
        if len(items) > preview
        else f'<p><a href="{filename}">{filename}</a></p>'
    )
    return (
        f'<details><summary>{summary}</summary>'
        f'{more}<ul>{rows}</ul></details>'
    )
