"""Features individually predictive of TRS in the training cohorts, for trs_prediction.py.

For each ssGSEA score and scFoundation dimension, the 5-fold cross-validated R² of
TRS on that feature alone, in the training samples only, against 9,999 permutations
of TRS; Benjamini–Hochberg across both feature sets, kept at adjusted p < 0.05.
Writes the stage's file in `trs_prediction_selected_features`.

Example:
    pixi run trs_prediction_feature_selection primary
"""

import argparse

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy.stats import false_discovery_control
from sklearn.linear_model import LinearRegression
from sklearn.metrics import make_scorer, r2_score
from sklearn.model_selection import KFold, cross_val_score

from common import cohort, config, resolve, stage_samples

TRAINING_COHORTS = ["EORTC", "Discovery_Stead"]
SEED = 11
N_SPLITS = 5
N_PERMUTATIONS = 9_999
ALPHA = 0.05


def cv_r2(x: np.ndarray, y: np.ndarray, cv: KFold) -> float:
    return cross_val_score(LinearRegression(), x, y, cv=cv, scoring=make_scorer(r2_score)).mean()


def permutation_cv_r2(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """A feature's cross-validated R² and its permutation p-value.

    Permutation j shuffles TRS with seed 11 + j, so every feature meets the
    same permutations.
    """
    cv = KFold(n_splits=N_SPLITS, shuffle=True, random_state=SEED)
    observed = cv_r2(x, y, cv)
    null = np.array([
        cv_r2(x, np.random.default_rng(SEED + j).permutation(y), cv)
        for j in range(N_PERMUTATIONS)
    ])
    return observed, ((null >= observed).sum() + 1) / (N_PERMUTATIONS + 1)


def selection(features: dict[str, pd.DataFrame], y: np.ndarray) -> pd.DataFrame:
    rows = []
    for feature_set, table in features.items():
        results = Parallel(n_jobs=-1)(
            delayed(permutation_cv_r2)(table[[column]].to_numpy(), y) for column in table.columns
        )
        rows += [
            {"feature_set": feature_set, "feature": column, "score": score, "p": p}
            for column, (score, p) in zip(table.columns, results)
        ]
    table = pd.DataFrame(rows)
    table["p_adjusted"] = false_discovery_control(table["p"], method="bh")
    table["kept"] = table["p_adjusted"] < ALPHA
    return table


def main(stage: str) -> None:
    settings = config()
    samples = stage_samples(TRAINING_COHORTS, stage)
    y = samples[cohort("bulk", TRAINING_COHORTS[0])["trs_column"]].to_numpy()
    features = {
        name: pd.read_csv(resolve(path), index_col=0).loc[samples.index]
        for name, path in settings["trs_prediction_features"].items()
    }
    print(f"{stage}: {len(samples)} training samples, "
          + ", ".join(f"{name} {table.shape[1]}" for name, table in features.items()))

    table = selection(features, y)
    output = resolve(settings["trs_prediction_selected_features"][stage])
    output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output, index=False)
    print(table[table["kept"]].groupby("feature_set").size().to_string())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("stage", choices=["primary", "recurrent"])
    main(parser.parse_args().stage)
