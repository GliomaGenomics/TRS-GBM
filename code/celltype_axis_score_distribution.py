"""Distribution of per-cell JBS axis scores across cell types.

Every cell of a single-cell dataset is placed on the JBS axis and one
distribution is drawn per cell type, ordered by mean score, so the cell types sit
along the axis in the order they occupy it.

Example:
    pixi run celltype_axis_score_distribution nomura_2025
"""

import argparse

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from common import (cell_axis_scores, cell_colours, cohort, resolve, save_panel,
                    set_paper_theme)

PANEL_WIDTH = 4
ROW_HEIGHT = 0.5


def scored_cells(name: str) -> pd.DataFrame:
    """Every labelled cell of a dataset, with its cell type and its axis score."""
    entry = cohort("single_cell", name)
    data = ad.read_h5ad(resolve(entry["path"]))
    cells = pd.DataFrame({
        "CellType": data.obs[entry["cell_type_column"]],
        "Axis Score": cell_axis_scores(data, entry["pseudocount"]),
    })
    return cells[cells["CellType"].notna()]


def mean_score_order(cells: pd.DataFrame) -> np.ndarray:
    means = cells.groupby("CellType", observed=True)["Axis Score"].mean()
    return means.sort_values().dropna().index.values


def distribution_panel(cells: pd.DataFrame, order: np.ndarray, title: str):
    colours = cell_colours()
    fig, ax = plt.subplots(
        figsize=(PANEL_WIDTH, ROW_HEIGHT * cells["CellType"].nunique())
    )
    sns.boxplot(
        data=cells,
        x="Axis Score",
        y="CellType",
        hue="CellType",
        palette=colours,
        order=order,
        hue_order=order,
        # The hue only carries the palette; without dodge=False each box and
        # violin is offset into its own hue slot and the two no longer align.
        dodge=False,
        fill=False,
        fliersize=1,
        width=0.5,
        zorder=0,
        ax=ax,
    )
    sns.violinplot(
        data=cells,
        x="Axis Score",
        y="CellType",
        hue="CellType",
        palette=colours,
        order=order,
        hue_order=order,
        saturation=1,
        alpha=0.3,
        linewidth=0,
        zorder=1,
        dodge=False,
        ax=ax,
    )
    # Seaborn 0.13 draws the violins in plotting order but hands out their hue
    # colours in order of appearance in the data, so each is recoloured to
    # match its box.
    for violin, cell_type in zip(ax.collections[-len(order):], order):
        violin.set_facecolor(colours[cell_type])
        violin.set_alpha(0.3)

    ax.set_xlabel("JBS Axis Score")
    ax.set_title(title)
    return fig


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "dataset", help="single-cell dataset in config/datasets.yaml, e.g. nomura_2025"
    )
    name = parser.parse_args().dataset

    set_paper_theme()
    cells = scored_cells(name)
    fig = distribution_panel(cells, mean_score_order(cells),
                             cohort("single_cell", name)["label"])
    save_panel(fig, cells.rename_axis("Cell").reset_index(),
               f"celltype_axis_score_distribution_{name}")
    plt.close(fig)


if __name__ == "__main__":
    main()
