"""Polarisation of JBS axis scores between primary and recurrent tumours.

Densities of the axis scores (log2(FPKM + 0.01) · loadings) of the primary and the
recurrent tumours of one bulk cohort's included pairs, with one Gaussian fitted to
the primaries and a two-component Gaussian mixture to the recurrences. The source
data holds the scores and the mixture parameters.

Example:
    pixi run axis_score_polarisation Discovery
"""

import argparse

import anndata as ad
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from scipy.stats import norm
from sklearn.mixture import GaussianMixture

from common import (PALETTE, cohort, fpkm_axis_scores, resolve, save_panel, selected,
                    set_paper_theme)

STAGES = ("Primary", "Recurrent")
STAGE_COLOURS = {"Primary": "#0197f6", "Recurrent": PALETTE["Recurrent"]}

STAGE_COMPONENTS = {"Primary": 1, "Recurrent": 2}
# Each component's curve is the stage colour moved along all three channels by
# this much, lower-mean component first.
COMPONENT_SHADES = {1: (-50,), 2: (10, -50)}

N_INIT = 20
REG_COVAR = 1e-6
RANDOM_STATE = 11


def axis_scores_by_stage(name: str) -> pd.DataFrame:
    """Both tumours of every included pair in a bulk cohort, scored on the JBS axis."""
    entry = cohort("bulk", name)
    obs = ad.read_h5ad(resolve(entry["path"]), backed="r").obs

    patients = obs[entry["patient_column"]].astype(object)
    included = (
        selected(obs, entry["filters"])
        & obs[entry["include_column"]].to_numpy()
        & (obs[entry["pair_purity_column"]] >= entry["min_purity"])
    )
    paired = patients.isin(patients[included])

    stages = {entry["primary_value"]: "Primary", entry["recurrent_value"]: "Recurrent"}
    return pd.DataFrame({
        "Patient": patients[paired],
        "Stage": obs[entry["stage_column"]].astype(object)[paired].map(stages),
        "Axis score": fpkm_axis_scores(entry).loc[obs.index[paired]],
    }).reset_index(drop=True)


def mixture(x: np.ndarray, n_components: int) -> dict:
    """Gaussian mixture, components ordered by their mean."""
    fit = GaussianMixture(
        n_components, n_init=N_INIT, reg_covar=REG_COVAR, random_state=RANDOM_STATE
    ).fit(x.reshape(-1, 1))
    order = np.argsort(fit.means_.ravel())
    return {
        "weights": fit.weights_[order],
        "means": fit.means_.ravel()[order],
        "sds": np.sqrt(fit.covariances_.ravel())[order],
    }


def hex_shift(colour: str, shift: int) -> str:
    """Move a hex colour equally along all three channels, to shade the two components apart."""
    channels = (int(colour[start:start + 2], 16) for start in (1, 3, 5))
    return "#" + "".join(f"{min(255, max(0, channel + shift)):02x}" for channel in channels)


def polarisation_panel(scores: pd.DataFrame):
    # Axes of 4.5 × 3 in, 8 pt tick labels and 12 pt text: the proportions of
    # the paper's panel.
    fig, ax = plt.subplots(figsize=(5.8, 3.9))
    # Key holds the patient for a score and the component (1 = lower mean) for
    # a mixture parameter.
    rows = []

    for stage in STAGES:
        at_stage = scores[scores["Stage"] == stage]
        x = at_stage["Axis score"].to_numpy()
        n_components = STAGE_COMPONENTS[stage]
        fit = mixture(x, n_components)
        colour = STAGE_COLOURS[stage]
        sns.kdeplot(x, bw_adjust=1, color=colour, fill=True, alpha=0.2,
                    linewidth=0, ax=ax, zorder=0)

        support = np.linspace(x.min() - 10, x.max() + 10, 1000)
        for component, shift in enumerate(COMPONENT_SHADES[n_components]):
            density = fit["weights"][component] * norm.pdf(
                support, fit["means"][component], fit["sds"][component]
            )
            ax.plot(support, density, linestyle="--",
                    color=hex_shift(colour, shift), zorder=1)

        rows += [{"Quantity": "Axis score", "Stage": stage, "Key": patient, "Value": score}
                 for patient, score in zip(at_stage["Patient"], x)]
        for parameter in ("weights", "means", "sds"):
            rows += [{"Quantity": f"GMM {parameter}", "Stage": stage,
                      "Key": component + 1, "Value": value}
                     for component, value in enumerate(fit[parameter])]

    handles = [
        mpatches.Patch(color=STAGE_COLOURS["Primary"], label="Primary GBM", alpha=0.2),
        mpatches.Patch(color=STAGE_COLOURS["Recurrent"], label="Recurrent GBM", alpha=0.2),
        Line2D([0], [0], color="grey", linestyle="--", label="Gaussian Mixture Model"),
    ]
    ax.legend(handles=handles, frameon=False, loc="upper left", fontsize=12,
              handlelength=1, handletextpad=0.4, labelspacing=0, borderaxespad=0)
    ax.set_xlabel("JBS Axis Score", fontsize=12)
    # Density units carry no meaning here: ticks are drawn, values and label are not.
    ax.set_ylabel("")
    ax.tick_params(axis="y", labelleft=False)

    return fig, pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("cohort", help="bulk cohort in config/datasets.yaml, e.g. Discovery")
    name = parser.parse_args().cohort

    set_paper_theme()
    fig, table = polarisation_panel(axis_scores_by_stage(name))
    save_panel(fig, table, f"axis_score_polarisation_{name}")
    plt.close(fig)


if __name__ == "__main__":
    main()
