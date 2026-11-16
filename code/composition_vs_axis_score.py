"""Sample composition along the JBS axis.

One stacked bar per sample, ordered by a pseudobulk JBS axis score drawn over the
bars, for every combination of

    compositions  population, microenvironment, cancer_cell, combined, axis_order
    scores        whole_tumour, microenvironment, cancer_cell

Groups are stacked and coloured by their mean per-cell axis score, so the cells are
read as well as the pseudobulks. Views drawn from the cancer cells alone keep
samples with purity above MIN_PURITY.

Writes each view to composition_vs_axis_score_<dataset>/<composition>_by_<score>_score
and all of them as one grid, composition_vs_axis_score_<dataset>.

Example:
    pixi run composition_vs_axis_score nomura_2025
"""

import argparse

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from common import cell_axis_scores, cohort, pseudobulk_table, resolve, save_panel, set_paper_theme

FIGSIZE = (6, 3.75)

# Each ramp as (colour map, start, end). The palest end of each is left off so
# that every segment stands out against the white page.
CANCER_CELL_RAMP = ("YlOrRd", 0.3, 1)
MICROENVIRONMENT_RAMP = ("YlGnBu", 0.3, 1)
AXIS_ORDER_RAMP = ("magma_r", 0.15, 0.9)
# The malignant population, in the dark red that ends the states' ramp.
MALIGNANT_POPULATION_COLOUR = "#800126"

# Fraction of a sample's cells that are malignant, at or below which the sample is
# too impure for its malignant state composition or cancer-cell score to be read.
MIN_PURITY = 0.1

SCORES = {
    "whole_tumour": ("All", "Whole-tumour JBS Axis Score"),
    "microenvironment": ("Non-malignant", "Microenvironment JBS Axis Score"),
    "cancer_cell": ("Malignant", "Cancer-cell JBS Axis Score"),
}


def axis_positions(entry: dict) -> dict[str, pd.Series]:
    """Mean per-cell axis score of each cell type, malignant state and cell population."""
    data = ad.read_h5ad(resolve(entry["path"]))
    cell_types = data.obs[entry["cell_type_column"]].astype("object")
    cells = pd.DataFrame({
        "cell_type": cell_types,
        "state": data.obs[entry["malignant_state_column"]].astype("object"),
        "population": cell_types.map(entry["cell_populations"]),
        "Axis Score": cell_axis_scores(data, entry["pseudocount"]),
    })
    return {group: cells.groupby(group, observed=True)["Axis Score"].mean()
            for group in ("cell_type", "state", "population")}


def position_colour_scale(positions: np.ndarray, ramp: tuple):
    """Colour for an axis position: the given positions are spaced evenly along the
    ramp in rank order, and any position between two of them is interpolated."""
    name, start, end = ramp
    ranked = np.sort(positions)
    steps = np.linspace(start, end, len(ranked))
    colour_map = plt.get_cmap(name)
    return lambda position: colour_map(np.interp(position, ranked, steps))


def composition_views(entry: dict, positions: dict[str, pd.Series]) -> dict:
    """Each composition: the pseudobulk its proportions come from (so that they sum
    to one), the columns it stacks from the bottom up with their labels and
    colours, its title, and legend headings to insert before the bars at the given
    positions."""
    states = positions["state"][entry["malignant_state_order"]].sort_values()
    cell_types = positions["cell_type"][entry["tme_order"]].sort_values()
    cancer_colour = position_colour_scale(states.to_numpy(), CANCER_CELL_RAMP)
    microenvironment_colour = position_colour_scale(cell_types.to_numpy(),
                                                    MICROENVIRONMENT_RAMP)
    axis_order_colour = position_colour_scale(
        np.concatenate([states, cell_types]), AXIS_ORDER_RAMP
    )

    def cancer(state: str) -> tuple:
        return cancer_colour(states[state])

    def microenvironment(cell_type: str) -> tuple:
        return microenvironment_colour(cell_types[cell_type])

    non_malignant = [name for name in entry["population_order"] if name != "Malignant"]
    populations = positions["population"][non_malignant].sort_values()

    pooled = pd.concat([
        pd.DataFrame({"column": "comp_malignant_" + states.index,
                      "label": states.index + " (malignant)", "position": states.values}),
        pd.DataFrame({"column": "comp_" + cell_types.index,
                      "label": cell_types.index, "position": cell_types.values}),
    ]).sort_values("position")

    return {
        "population": {
            "population": "All", "title": "Population", "headings": {},
            "parts": {
                "comp_pop_Malignant": ("Malignant", MALIGNANT_POPULATION_COLOUR),
                **{"comp_pop_" + name: (name, microenvironment_colour(position))
                   for name, position in populations.items()},
            },
        },
        "microenvironment": {
            "population": "Non-malignant", "title": "Microenvironment", "headings": {},
            "parts": {"comp_" + name: (name, microenvironment(name))
                      for name in cell_types.index},
        },
        "cancer_cell": {
            "population": "Malignant", "title": "Cancer-cell", "headings": {},
            "parts": {"malignant_comp_" + name: (name, cancer(name))
                      for name in states.index},
        },
        "combined": {
            "population": "All", "title": "Cancer-cell and microenvironment",
            "headings": {0: "Cancer cells", len(states): "Microenvironment"},
            "parts": {
                **{"comp_malignant_" + name: (name, cancer(name)) for name in states.index},
                **{"comp_" + name: (name, microenvironment(name)) for name in cell_types.index},
            },
        },
        "axis_order": {
            "population": "All", "title": "Axis-ordered", "headings": {},
            "parts": {row.column: (row.label, axis_order_colour(row.position))
                      for row in pooled.itertuples()},
        },
    }


def composition_table(obs: pd.DataFrame, entry: dict, view: dict,
                      score: str) -> pd.DataFrame:
    """The proportions of one composition, one row per sample, in increasing order
    of one axis score, which is the last column."""
    scored_by, score_label = SCORES[score]

    if "Malignant" in (view["population"], scored_by):
        whole = obs[obs["Population"] == "All"]
        pure = whole.loc[whole["Purity"] > MIN_PURITY, "Sample"]
        obs = obs[obs["Sample"].isin(pure)]

    bars = obs[obs["Population"] == view["population"]].set_index("Sample")
    scores = obs[obs["Population"] == scored_by].set_index("Sample")["Axis Score"]
    return (bars[list(view["parts"])].assign(**{score_label: scores})
            .sort_values(score_label))


def stacked_bars(ax, table: pd.DataFrame, parts: dict):
    """Stacked composition bars with the axis score on a second y axis, which is
    returned.

    Sample names are left off the x axis; at this many samples they are unreadable
    and the source data table carries them.
    """
    *columns, score_label = table.columns

    bottoms = np.zeros(len(table))
    for column in columns:
        values = table[column].to_numpy()
        label, colour = parts[column]
        ax.bar(table.index, values, bottom=bottoms, label=label,
               color=colour, edgecolor="none", width=1.0)
        bottoms += values

    ax.margins(x=0)
    ax.set_xticks([])
    ax.set_xlabel("Sample")
    ax.set_ylabel("Proportion")
    ax.set_ylim(0, 1)
    ax.tick_params(axis="y", labelsize=9)

    ax2 = ax.twinx()
    ax2.plot(table.index, table[score_label], color="black", marker="none",
             linestyle="-", linewidth=1, label="JBS Axis Score")
    ax2.set_ylabel(score_label, color="black")
    ax2.tick_params(axis="y", labelcolor="black", labelsize=9)
    ax2.spines["right"].set_visible(True)
    return ax2


def draw_legend(ax, ax2, headings: dict, **placement) -> None:
    """One legend for the bars, bottom first, and the score line, with a bold
    heading before the bar at each position in headings."""
    bars, bar_labels = ax.get_legend_handles_labels()
    lines, line_labels = ax2.get_legend_handles_labels()
    handles, labels = [], []
    for position, (bar, label) in enumerate(zip(bars, bar_labels)):
        if position in headings:
            handles.append(Line2D([], [], linestyle="none"))
            labels.append(headings[position])
        handles.append(bar)
        labels.append(label)

    legend = ax2.legend(handles + lines, labels + line_labels, **placement)
    for text in legend.get_texts():
        if text.get_text() in headings.values():
            text.set_fontweight("bold")


def composition_panel(table: pd.DataFrame, view: dict, title: str):
    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax2 = stacked_bars(ax, table, view["parts"])
    ax.set_title(title)
    if len(view["parts"]) > 10:
        # Laid out before the legend is added, so that the legend hangs beside
        # the axes rather than shrinking them.
        fig.tight_layout()
        draw_legend(ax, ax2, view["headings"], loc="upper left",
                    bbox_to_anchor=(1.2, 1), frameon=False)
    else:
        draw_legend(ax, ax2, view["headings"], loc="upper left",
                    ncols=int((len(table.columns) + 1) / 4), framealpha=0.75)
        fig.tight_layout()
    return fig


def composition_grid(tables: dict, views: dict, title: str):
    """Every composition against every score: compositions in rows, scores in
    columns, one legend per row beside its last panel.

    The panels are laid out before the legends are added, so that a legend taller
    than its row hangs beside the panels rather than pushing the rows apart.
    """
    fig, axes = plt.subplots(len(views), len(SCORES), squeeze=False,
                             figsize=(FIGSIZE[0] * len(SCORES), FIGSIZE[1] * len(views)))
    row_ends = []
    for row, (composition, view) in enumerate(views.items()):
        for column, score in enumerate(SCORES):
            ax = axes[row, column]
            table = tables[composition, score]
            ax2 = stacked_bars(ax, table, view["parts"])
            ax.set_title(f"{view['title']} composition (n = {len(table)})")
        row_ends.append((ax, ax2, view["headings"]))

    fig.tight_layout()
    for ax, ax2, headings in row_ends:
        draw_legend(ax, ax2, headings, loc="upper left", bbox_to_anchor=(1.2, 1.05),
                    frameon=False, fontsize=8, labelspacing=0.3)
    fig.suptitle(title, y=1, va="bottom")
    return fig


def grid_source_data(tables: dict) -> pd.DataFrame:
    """Every panel's values in one long table, one row per sample and bar segment."""
    long = []
    for (composition, score), table in tables.items():
        *parts, score_label = table.columns
        melted = table.reset_index().melt(
            id_vars=["Sample", score_label], value_vars=parts,
            var_name="Part", value_name="Proportion",
        )
        long.append(melted.rename(columns={score_label: "Axis Score"})
                    .assign(Composition=composition, Score=score))
    columns = ["Composition", "Score", "Sample", "Axis Score", "Part", "Proportion"]
    return pd.concat(long)[columns]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "dataset", help="single-cell dataset in config/datasets.yaml, e.g. nomura_2025"
    )
    arguments = parser.parse_args()

    set_paper_theme()
    entry = cohort("single_cell", arguments.dataset)
    obs = pseudobulk_table(arguments.dataset)
    views = composition_views(entry, axis_positions(entry))
    name = f"composition_vs_axis_score_{arguments.dataset}"

    tables = {}
    for composition, view in views.items():
        for score in SCORES:
            table = composition_table(obs, entry, view, score)
            fig = composition_panel(table, view, entry["label"])
            save_panel(fig, table.reset_index(),
                       f"{name}/{composition}_by_{score}_score")
            plt.close(fig)
            tables[composition, score] = table

    fig = composition_grid(tables, views, entry["label"])
    save_panel(fig, grid_source_data(tables), name)
    plt.close(fig)


if __name__ == "__main__":
    main()
