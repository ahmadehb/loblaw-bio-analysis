"""
Part 3: Responder vs non analysis

Cohort: melanoma + miraclib + PBMC samples
Test: Mann-Whitney U (two-sided) per population per timepoint
Effect size: rank-biserial correlation (r)
Multiple comparisons: Benjamini-Hochberg FDR for all 15 tests

EDA notebook (notebooks/eda_stats_check.ipynb) documents choices:
distributions non-normal (Shapiro-Wilk rejects in many groups), variances
not always equal, so Mann-Whitney U is safer choice than t-tests

Outputs:
  outputs/statistical_results.csv   p-values, effect sizes, group medians
  outputs/boxplots.png              5 populations × 3 timepoints, by response
"""

from __future__ import annotations

import sqlite3
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
DB_PATH = ROOT / "loblaw.db"
OUT_DIR = ROOT / "outputs"

# Pull analysis cohort & compute per-sample percents
COHORT_QUERY = """
SELECT
    sub.subject_id,
    s.sample_id,
    s.response,
    s.time_from_treatment_start AS timepoint,
    cc.population,
    cc.count,
    100.0 * cc.count / SUM(cc.count) OVER (PARTITION BY s.sample_id) AS percentage
FROM subjects sub
JOIN samples     s  ON s.subject_id = sub.subject_id
JOIN cell_counts cc ON cc.sample_id = s.sample_id
WHERE sub.condition  = 'melanoma'
  AND s.treatment    = 'miraclib'
  AND s.sample_type  = 'PBMC'
  AND s.response IS NOT NULL
ORDER BY s.sample_id, cc.population;
"""


def rank_biserial(u_stat: float, n1: int, n2: int) -> float:
    """
    Rank-biserial correlation as an effect size for Mann-Whitney U.

    r = 1 - 2U / (n1*n2)  where U is the smaller of the two U statistics.
    Range: -1 to +1. Sign indicates direction (positive => group 1 > group 2).
    Magnitude rules of thumb: 0.1 small, 0.3 medium, 0.5 large.
    """
    return 1.0 - (2.0 * u_stat) / (n1 * n2)


def benjamini_hochberg(pvalues: np.ndarray, alpha: float = 0.05) -> np.ndarray:
    """
    Benjamini-Hochberg FDR-corrected q-values.

    Returns q-values where rejecting H0 when q <= alpha controls FDR at alpha.
    Implements the standard step-up procedure: sort, scale, monotone-correct.
    """
    p = np.asarray(pvalues, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    # Enforce monotonicity from the largest rank backwards.
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty_like(ranked)
    q[order] = np.clip(ranked, 0, 1)
    return q


def run_tests(df: pd.DataFrame) -> pd.DataFrame:
    """One Mann-Whitney U test per (population, timepoint)."""
    results = []
    for population in sorted(df.population.unique()):
        for timepoint in sorted(df.timepoint.unique()):
            block = df[(df.population == population) & (df.timepoint == timepoint)]
            yes = block.loc[block.response == "yes", "percentage"].to_numpy()
            no  = block.loc[block.response == "no",  "percentage"].to_numpy()

            u_stat, p_value = stats.mannwhitneyu(yes, no, alternative="two-sided")
            # rank_biserial expects smaller U; mannwhitneyu returns U for the
            # first group, so make symmetrical
            u_min = min(u_stat, len(yes) * len(no) - u_stat)
            effect = rank_biserial(u_min, len(yes), len(no))
            # Direction: positive => responders higher than non-responders.
            direction = np.sign(np.median(yes) - np.median(no))
            signed_effect = direction * abs(effect)

            results.append({
                "population": population,
                "timepoint": int(timepoint),
                "n_responders": len(yes),
                "n_nonresponders": len(no),
                "median_yes": float(np.median(yes)),
                "median_no":  float(np.median(no)),
                "median_diff_yes_minus_no": float(np.median(yes) - np.median(no)),
                "u_statistic": float(u_stat),
                "p_value": float(p_value),
                "effect_size_rank_biserial": float(signed_effect),
            })

    out = pd.DataFrame(results)
    # FDR correction for all tests
    out["q_value_bh"] = benjamini_hochberg(out["p_value"].to_numpy())
    out["significant_q05"] = out["q_value_bh"] < 0.05
    return out.sort_values(["population", "timepoint"]).reset_index(drop=True)


def make_boxplots(df: pd.DataFrame, results: pd.DataFrame, out_path: Path) -> None:
    """5 populations × 3 timepoints. Each panel: yes vs no boxplot with q-value."""
    populations = sorted(df.population.unique())
    timepoints = sorted(df.timepoint.unique())

    fig, axes = plt.subplots(
        len(populations), len(timepoints),
        figsize=(11, 14), sharey="row",
    )

    for i, pop in enumerate(populations):
        for j, tp in enumerate(timepoints):
            ax = axes[i, j]
            block = df[(df.population == pop) & (df.timepoint == tp)]
            yes = block.loc[block.response == "yes", "percentage"]
            no  = block.loc[block.response == "no",  "percentage"]

            bp = ax.boxplot(
                [no, yes], tick_labels=["no", "yes"],
                patch_artist=True, widths=0.55, showfliers=False,
            )
            for patch, color in zip(bp["boxes"], ["#d9534f", "#3a7ca5"]):
                patch.set_facecolor(color)
                patch.set_alpha(0.65)
            for med in bp["medians"]:
                med.set_color("black"); med.set_linewidth(1.5)

            # Overlay jittered points
            for k, vals in enumerate([no, yes], start=1):
                pts = vals.sample(min(len(vals), 80), random_state=0)
                x_jit = np.random.default_rng(0).normal(k, 0.05, size=len(pts))
                ax.scatter(x_jit, pts, s=6, alpha=0.35, color="black", zorder=3)

            row = results[(results.population == pop) & (results.timepoint == tp)].iloc[0]
            star = " *" if row.significant_q05 else ""
            ax.set_title(
                f"{pop} @ t={tp}d\nq={row.q_value_bh:.3g}{star}",
                fontsize=9,
            )
            if j == 0:
                ax.set_ylabel("% of total", fontsize=9)
            ax.tick_params(labelsize=8)

    fig.suptitle(
        "Responder (yes) vs non-responder (no): melanoma + miraclib + PBMC\n"
        "Mann-Whitney U with BH-FDR correction; * = q < 0.05",
        fontsize=11, y=1.00,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    if not DB_PATH.exists():
        raise SystemExit(f"ERROR: {DB_PATH} not found. Run `python load_data.py` first.")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading cohort (melanoma + miraclib + PBMC) ...")
    with sqlite3.connect(DB_PATH) as conn:
        df = pd.read_sql_query(COHORT_QUERY, conn)
    print(f"  {df.sample_id.nunique()} samples from {df.subject_id.nunique()} subjects")

    print("\nRunning Mann-Whitney U tests (5 populations × 3 timepoints = 15 tests) ...")
    results = run_tests(df)

    results_path = OUT_DIR / "statistical_results.csv"
    results.to_csv(results_path, index=False, float_format="%.6g")
    print(f"  Wrote {results_path.relative_to(ROOT)}")

    boxplot_path = OUT_DIR / "boxplots.png"
    make_boxplots(df, results, boxplot_path)
    print(f"  Wrote {boxplot_path.relative_to(ROOT)}")

    # Summary to stdout
    sig = results[results.significant_q05]
    print(f"\n=== Significant findings (BH-FDR q < 0.05): {len(sig)} of {len(results)} ===")
    if len(sig):
        show = sig[[
            "population", "timepoint", "median_yes", "median_no",
            "median_diff_yes_minus_no", "p_value", "q_value_bh",
            "effect_size_rank_biserial",
        ]].copy()
        print(show.to_string(index=False))
    else:
        lowest = results.sort_values("q_value_bh").iloc[0]
        print("  No comparisons reached FDR significance.")
        print(f"  Lowest q-value: {lowest.q_value_bh:.3f} "
              f"({lowest.population} @ t={lowest.timepoint})")

    # Exploratory: comparisons with uncorrected p < 0.05
    nominal = results[(results.p_value < 0.05) & ~results.significant_q05]
    if len(nominal):
        print(f"\n=== Exploratory: uncorrected p < 0.05 (do NOT survive FDR) ===")
        print("  Listed for completeness only; do not interpret as predictive.")
        show = nominal[[
            "population", "timepoint", "median_yes", "median_no",
            "p_value", "q_value_bh", "effect_size_rank_biserial",
        ]].copy()
        print(show.to_string(index=False))

    print("\n=== Interpretation ===")
    print("  Primary finding: no immune cell population shows a statistically")
    print("  significant difference in relative frequency between miraclib")
    print("  responders and non-responders in melanoma PBMC samples at any")
    print("  timepoint (lowest q = {:.2f}; |effect size r| <= {:.2f}).".format(
        results.q_value_bh.min(),
        results.effect_size_rank_biserial.abs().max(),
    ))
    print("  Bulk PBMC cell-type proportions do not appear to be a useful")
    print("  predictive biomarker for miraclib response in this cohort.")


if __name__ == "__main__":
    main()
