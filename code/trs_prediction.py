"""TRS predicted from one stage's tumour expression, tested blind on held-out centres.

Trained on one stage of the included pairs of EORTC and Discovery_Stead, tested on
Discovery_Heldout, with the ssGSEA scores and scFoundation dimensions that
trs_prediction_feature_selection.py kept. Models: ElasticNet, and a random forest
on the k features of highest mutual information with TRS, tuned by BayesSearchCV
(60 iterations, R², 5 × 5-fold cross-validation); nested cross-validation for
training performance; seed 11.

Writes the panel trs_prediction_from_<stage> and the performance table
figures/tables/trs_prediction_from_<stage>.csv.

Example:
    pixi run trs_prediction primary
"""

import argparse
from functools import partial

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import yaml
from matplotlib.legend_handler import HandlerTuple
from matplotlib.lines import Line2D
from scipy.stats import sem
from sklearn.ensemble import RandomForestRegressor
from sklearn.feature_selection import SelectKBest, mutual_info_regression
from sklearn.linear_model import ElasticNet
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import KFold, RepeatedKFold, cross_validate
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from skopt import BayesSearchCV
from skopt.space import Integer, Real

from common import FIGURES, cohort, config, resolve, save_panel, stage_samples

TRAINING_COHORTS = ["EORTC", "Discovery_Stead"]
TEST_COHORT = "Discovery_Heldout"
FEATURE_SETS = ["ssGSEA", "scFoundation"]
MODELS = ["ElasticNet", "RandomForest"]
SEED = 11
N_ITERATIONS = 60
TABLES = FIGURES / "tables"


def set_theme() -> None:
    sns.set_theme(
        style="ticks",
        context="paper",
        font_scale=1,
        font="sans-serif",
        rc={
            "axes.linewidth": 1.2,
            "xtick.major.width": 1.2,
            "ytick.major.width": 1.2,
            "savefig.transparent": True,
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
        },
    )


def centre_colours() -> dict:
    with open(resolve("config/palette.yaml")) as handle:
        return yaml.safe_load(handle)["centre"]


def kept_features(stage: str, feature_set: str) -> list[str]:
    selection = pd.read_csv(resolve(config()["trs_prediction_selected_features"][stage]))
    kept = selection[selection["kept"] & (selection["feature_set"] == feature_set)]
    return kept["feature"].tolist()


def search(model: str, n_features: int, inner_cv: RepeatedKFold) -> BayesSearchCV:
    """The hyperparameter search for one model, ready to fit."""
    if model == "ElasticNet":
        pipeline = Pipeline([
            ("scaler", StandardScaler()),
            ("model", ElasticNet(max_iter=1_000_000, random_state=SEED)),
        ])
        space = {
            "model__alpha": Real(1e-4, 1e2, prior="log-uniform"),
            "model__l1_ratio": Real(0, 1, prior="uniform"),
        }
    else:
        pipeline = Pipeline([
            # Mutual information is estimated with added noise; seeded so runs repeat.
            ("selector", SelectKBest(score_func=partial(mutual_info_regression, random_state=SEED))),
            ("model", RandomForestRegressor(criterion="squared_error", random_state=SEED)),
        ])
        space = {
            "model__max_depth": Integer(1, 60),
            "model__n_estimators": Integer(10, 300),
            "model__min_samples_leaf": Integer(1, 25),
            "model__min_samples_split": Integer(2, 50),
            "model__max_leaf_nodes": Integer(2, 50),
            "selector__k": Integer(1, n_features),
        }
    return BayesSearchCV(
        estimator=pipeline,
        search_spaces=space,
        n_iter=N_ITERATIONS,
        cv=inner_cv,
        scoring="r2",
        n_jobs=-1,
        random_state=SEED,
    )


def evaluate(model: str, X_train: np.ndarray, y_train: np.ndarray, X_test: np.ndarray,
             y_test: np.ndarray) -> tuple[dict, np.ndarray, np.ndarray]:
    """Nested cross-validation, then the refitted model's training and test predictions."""
    outer_cv = KFold(n_splits=5, shuffle=True, random_state=SEED)
    inner_cv = RepeatedKFold(n_splits=5, n_repeats=5, random_state=SEED)
    tuned = search(model, X_train.shape[1], inner_cv)

    nested = cross_validate(
        tuned, X=X_train, y=y_train, cv=outer_cv,
        scoring={"RMSE": "neg_root_mean_squared_error", "R2": "r2"},
    )
    tuned.fit(X_train, y_train)
    train_predicted = tuned.predict(X_train)
    test_predicted = tuned.predict(X_test)

    performance = {
        "Nested CV RMSE": -nested["test_RMSE"].mean(),
        "Nested CV RMSE SEM": sem(-nested["test_RMSE"]),
        "Nested CV R²": nested["test_R2"].mean(),
        "Nested CV R² SEM": sem(nested["test_R2"]),
        "Training R²": r2_score(y_train, train_predicted),
        "Training RMSE": np.sqrt(mean_squared_error(y_train, train_predicted)),
        "Test R²": r2_score(y_test, test_predicted),
        "Test RMSE": np.sqrt(mean_squared_error(y_test, test_predicted)),
        "Features kept": (tuned.best_estimator_.named_steps["selector"].k
                          if model == "RandomForest" else None),
    }
    return performance, train_predicted, test_predicted


def draw(ax, run: pd.DataFrame, colours: dict, title: str) -> None:
    train = run[run["Set"] == "Train"]
    test = run[run["Set"] == "Test"]
    sns.scatterplot(
        y=train["TRS"].to_numpy(), x=train["Predicted TRS"].to_numpy(), color="lightgrey",
        s=20, label="Train", linewidth=0.5, edgecolor="darkgrey", alpha=1, ax=ax,
    )
    sns.scatterplot(
        y=test["TRS"].to_numpy(), x=test["Predicted TRS"].to_numpy(),
        hue=test["Centre"].to_numpy(), palette=colours, s=20, linewidth=0.5,
        edgecolor="black", alpha=1, legend=False, ax=ax,
    )
    limits = [np.min([ax.get_xlim(), ax.get_ylim()]), np.max([ax.get_xlim(), ax.get_ylim()])]
    ax.plot(limits, limits, linestyle="--", color="black", alpha=1)
    ax.set_ylabel("True TRS")
    ax.set_xlabel("Predicted TRS")
    ax.set_title(title)

    train_proxy = Line2D([], [], marker="o", color="darkgrey", markerfacecolor="lightgrey",
                         linestyle="None", markersize=5, alpha=1, label="Train")
    test_proxies = tuple(
        Line2D([], [], marker="o", color=colours[centre], markerfacecolor=colours[centre],
               markeredgecolor="black", linestyle="None", markersize=5, alpha=1)
        for centre in colours if centre in set(test["Centre"])
    )
    ax.legend(handles=[train_proxy, test_proxies], labels=["Train", "Test"],
              handler_map={tuple: HandlerTuple(ndivide=None)}, loc="upper left")


def main(stage: str) -> None:
    settings = config()
    entry = cohort("bulk", TEST_COHORT)
    trs, centre = entry["trs_column"], entry["centre_column"]
    train = stage_samples(TRAINING_COHORTS, stage)
    test = stage_samples([TEST_COHORT], stage)
    y_train, y_test = train[trs].to_numpy(), test[trs].to_numpy()

    runs, performance = [], []
    for feature_set in FEATURE_SETS:
        features = pd.read_csv(resolve(settings["trs_prediction_features"][feature_set]), index_col=0)
        features = features.loc[:, kept_features(stage, feature_set)]
        X_train = features.loc[train.index].to_numpy()
        X_test = features.loc[test.index].to_numpy()
        print(f"{stage} {feature_set}: train {X_train.shape}, test {X_test.shape}")

        for model in MODELS:
            scores, train_predicted, test_predicted = evaluate(model, X_train, y_train, X_test, y_test)
            performance.append({"Features": feature_set, "Model": model, **scores})
            for samples, y, predicted, label in [(train, y_train, train_predicted, "Train"),
                                                 (test, y_test, test_predicted, "Test")]:
                runs.append(pd.DataFrame({
                    "Sample": samples.index,
                    "Patient": samples[entry["patient_column"]].to_numpy(),
                    "Centre": samples[centre].astype(str).to_numpy(),
                    "Set": label,
                    "Features": feature_set,
                    "Model": model,
                    "TRS": y,
                    "Predicted TRS": predicted,
                }))
    table = pd.concat(runs, ignore_index=True)
    performance = pd.DataFrame(performance)

    set_theme()
    colours = centre_colours()
    fig, axes = plt.subplots(len(MODELS), len(FEATURE_SETS), figsize=(6, 6))
    for row, model in enumerate(MODELS):
        for column, feature_set in enumerate(FEATURE_SETS):
            run = table[(table["Model"] == model) & (table["Features"] == feature_set)]
            draw(axes[row, column], run, colours, f"{feature_set} {model}")
    fig.tight_layout()

    name = f"trs_prediction_from_{stage}"
    save_panel(fig, table, name)
    TABLES.mkdir(parents=True, exist_ok=True)
    performance.to_csv(TABLES / f"{name}.csv", index=False)
    print(performance.to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("stage", choices=["primary", "recurrent"])
    main(parser.parse_args().stage)
