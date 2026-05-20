"""
Streamlit dashboard for the Loblaw Bio cell-count analysis

Layout:
  Sidebar  — global filters (condition, treatment, sample type, timepoint) that
             apply to Overview, Frequencies, and Baseline-subset tabs.
  Tabs     — 1. Overview        cohort sizes and demographics
             2. Frequencies     Part 2 frequency table, faceted by population
             3. Response stats  Part 3 (always melanoma + miraclib + PBMC,
                                the required comparison) — boxplots, stats,
                                interpretation. Optional "explore other
                                cohorts" toggle re-runs the same analysis on
                                a user-defined subset for exploration.
             4. Baseline subset Part 4 baseline counts and the male-responder
                                B-cell average.

Run from the repo root with:  streamlit run dashboard/app.py
"""

from __future__ import annotations

import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
DB_PATH = ROOT / "loblaw.db"

st.set_page_config(
    page_title="Loblaw Bio — Cell Count Analysis",
    page_icon="🧬",
    layout="wide",
)


# Data access

@st.cache_data
def load_long() -> pd.DataFrame:
    """Load the full long-format frame once and cache it.

    Returns one row per (sample, population) with denormalized subject and
    sample metadata, plus a per-sample percentage computed in SQL.
    """
    if not DB_PATH.exists():
        st.error(
            f"Database not found at `{DB_PATH}`.\n\n"
            "Run `python load_data.py` from the repo root before launching the dashboard."
        )
        st.stop()

    with sqlite3.connect(DB_PATH) as conn:
        df = pd.read_sql_query(
            """
            SELECT
                sub.subject_id, sub.project_id, sub.condition, sub.age, sub.sex,
                s.sample_id, s.sample_type, s.treatment, s.response,
                s.time_from_treatment_start AS timepoint,
                cc.population, cc.count,
                100.0 * cc.count /
                    SUM(cc.count) OVER (PARTITION BY s.sample_id) AS percentage
            FROM subjects sub
            JOIN samples s ON s.subject_id = sub.subject_id
            JOIN cell_counts cc ON cc.sample_id = s.sample_id
            """,
            conn,
        )
    return df



# Statistical helpers


def benjamini_hochberg(pvalues: np.ndarray) -> np.ndarray:
    """BH step-up FDR; returns q-values (same shape as input)."""
    p = np.asarray(pvalues, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty_like(ranked)
    q[order] = np.clip(ranked, 0, 1)
    return q


def run_response_tests(df: pd.DataFrame) -> pd.DataFrame:
    """Mann-Whitney U per (population, timepoint) on the given filtered frame.

    Expects df to already be limited to a single condition + treatment + sample
    type, with response in {'yes','no'} and timepoint, population, percentage
    columns.
    """
    rows = []
    for pop in sorted(df.population.unique()):
        for tp in sorted(df.timepoint.unique()):
            block = df[(df.population == pop) & (df.timepoint == tp)]
            yes = block.loc[block.response == "yes", "percentage"].to_numpy()
            no = block.loc[block.response == "no", "percentage"].to_numpy()
            if len(yes) < 3 or len(no) < 3:
                # Not enough data to run meaningful test
                rows.append({
                    "population": pop, "timepoint": int(tp),
                    "n_yes": len(yes), "n_no": len(no),
                    "median_yes": np.nan, "median_no": np.nan,
                    "p_value": np.nan,
                })
                continue
            u, p = stats.mannwhitneyu(yes, no, alternative="two-sided")
            rows.append({
                "population": pop, "timepoint": int(tp),
                "n_yes": len(yes), "n_no": len(no),
                "median_yes": float(np.median(yes)),
                "median_no": float(np.median(no)),
                "p_value": float(p),
            })
    out = pd.DataFrame(rows)
    valid = out.p_value.notna()
    out["q_value_bh"] = np.nan
    if valid.any():
        out.loc[valid, "q_value_bh"] = benjamini_hochberg(out.loc[valid, "p_value"].to_numpy())
    out["significant_q05"] = out["q_value_bh"] < 0.05
    return out



# Sidebar filters


data = load_long()

with st.sidebar:
    st.title("🧬 Loblaw Bio")
    st.caption("Cell count analysis dashboard")
    st.divider()
    st.subheader("Global filters")
    st.caption("Apply to Overview, Frequencies, and Baseline tabs.")

    def multi(label: str, options: list, default: list | None = None):
        return st.multiselect(label, options, default if default is not None else options)

    conditions = sorted(data.condition.unique())
    treatments = sorted(data.treatment.unique())
    sample_types = sorted(data.sample_type.unique())
    timepoints = sorted(data.timepoint.unique())

    sel_condition = multi("Condition", conditions)
    sel_treatment = multi("Treatment", treatments)
    sel_sample_type = multi("Sample type", sample_types)
    sel_timepoint = multi("Timepoint (days)", timepoints)

    st.divider()
    st.caption(
        "**Note:** the Response Stats tab always runs the required "
        "melanoma + miraclib + PBMC comparison (Part 3 spec). "
        "Use the toggle on that tab to explore other cohorts."
    )

# Apply sidebar filters once; tabs read from `filtered`
filtered = data[
    data.condition.isin(sel_condition)
    & data.treatment.isin(sel_treatment)
    & data.sample_type.isin(sel_sample_type)
    & data.timepoint.isin(sel_timepoint)
]

# Tabs


tab_overview, tab_freq, tab_stats, tab_subset = st.tabs([
    "📊 Overview", "📈 Frequencies", "🔬 Response stats", "🎯 Baseline subset",
])


# Tab 1: Overview

with tab_overview:
    st.header("Cohort overview")

    if filtered.empty:
        st.warning("No samples match the current filters. Widen the selection in the sidebar.")
    else:
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Projects", filtered.project_id.nunique())
        c2.metric("Subjects", filtered.subject_id.nunique())
        c3.metric("Samples", filtered.sample_id.nunique())
        c4.metric("Cell counts", f"{len(filtered):,}")

        st.divider()

        col_left, col_right = st.columns(2)

        with col_left:
            st.subheader("Subjects by condition")
            sub_demo = filtered.drop_duplicates("subject_id")
            cond_counts = sub_demo.condition.value_counts().reset_index()
            cond_counts.columns = ["condition", "subjects"]
            fig = px.bar(
                cond_counts, x="condition", y="subjects",
                color="condition", text="subjects",
            )
            fig.update_layout(showlegend=False, height=320, margin=dict(t=10, b=10))
            st.plotly_chart(fig, use_container_width=True)

            st.subheader("Subjects by sex")
            sex_counts = sub_demo.sex.value_counts().reset_index()
            sex_counts.columns = ["sex", "subjects"]
            fig = px.bar(
                sex_counts, x="sex", y="subjects",
                color="sex", text="subjects",
                color_discrete_map={"M": "#3a7ca5", "F": "#d96a8a"},
            )
            fig.update_layout(showlegend=False, height=280, margin=dict(t=10, b=10))
            st.plotly_chart(fig, use_container_width=True)

        with col_right:
            st.subheader("Subjects by treatment × response")
            crosstab = (
                sub_demo.groupby(["treatment", "response"], dropna=False)
                .subject_id.nunique().reset_index(name="subjects")
            )
            crosstab["response"] = crosstab["response"].fillna("n/a")
            fig = px.bar(
                crosstab, x="treatment", y="subjects",
                color="response", text="subjects", barmode="group",
                color_discrete_map={"yes": "#3a7ca5", "no": "#d9534f", "n/a": "#888"},
            )
            fig.update_layout(height=320, margin=dict(t=10, b=10))
            st.plotly_chart(fig, use_container_width=True)

            st.subheader("Samples per project")
            proj_counts = (
                filtered.drop_duplicates("sample_id")
                .project_id.value_counts().reset_index()
            )
            proj_counts.columns = ["project_id", "samples"]
            fig = px.bar(
                proj_counts.sort_values("project_id"),
                x="project_id", y="samples", text="samples",
            )
            fig.update_layout(showlegend=False, height=280, margin=dict(t=10, b=10))
            st.plotly_chart(fig, use_container_width=True)


# Tab 2: Frequencies (Part 2)
with tab_freq:
    st.header("Per-sample cell population frequencies")
    st.caption(
        "Part 2 deliverable. Each row = one (sample, population) pair with the "
        "sample's total count and that population's relative frequency."
    )

    if filtered.empty:
        st.warning("No samples match the current filters.")
    else:
        freq_table = filtered[[
            "sample_id", "population", "count", "percentage",
        ]].rename(columns={"sample_id": "sample"})
        freq_table["total_count"] = filtered.groupby("sample_id")["count"].transform("sum")
        freq_table = freq_table[["sample", "total_count", "population", "count", "percentage"]]
        freq_table["percentage"] = freq_table["percentage"].round(4)
        freq_table = freq_table.sort_values(["sample", "population"]).reset_index(drop=True)

        c1, c2 = st.columns([1, 3])
        c1.metric("Rows", f"{len(freq_table):,}")
        c1.metric("Samples", f"{freq_table['sample'].nunique():,}")

        with c2:
            st.subheader("Distribution of relative frequencies by population")
            fig = px.box(
                filtered, x="population", y="percentage", color="population",
                points=False,
                category_orders={"population": sorted(filtered.population.unique())},
            )
            fig.update_layout(showlegend=False, height=360, margin=dict(t=10, b=10),
                              yaxis_title="% of total")
            st.plotly_chart(fig, use_container_width=True)

        st.subheader("Frequency table")
        st.dataframe(freq_table, use_container_width=True, height=400)

        st.download_button(
            "Download as CSV",
            freq_table.to_csv(index=False).encode("utf-8"),
            file_name="frequency_table_filtered.csv",
            mime="text/csv",
        )


# Tab 3: Response stats (Part 3)
with tab_stats:
    st.header("Responder vs non-responder analysis")
    st.caption(
        "Part 3 deliverable. By default, this tab always runs the required "
        "comparison: melanoma + miraclib + PBMC, responders (yes) vs "
        "non-responders (no), stratified by timepoint."
    )

    explore = st.toggle(
        "Explore other cohorts (overrides the default)",
        value=False,
        help="Re-run the same Mann-Whitney U + BH-FDR analysis on a different "
             "cohort. The default remains visible as the required deliverable.",
    )

    if explore:
        c1, c2, c3 = st.columns(3)
        ex_condition = c1.selectbox("Condition", sorted(data.condition.unique()), index=sorted(data.condition.unique()).index("melanoma"))
        ex_treatment = c2.selectbox("Treatment",
                                     [t for t in sorted(data.treatment.unique()) if t != "none"],
                                     index=0)
        ex_sample = c3.selectbox("Sample type", sorted(data.sample_type.unique()), index=sorted(data.sample_type.unique()).index("PBMC"))
        cohort = data[
            (data.condition == ex_condition)
            & (data.treatment == ex_treatment)
            & (data.sample_type == ex_sample)
            & data.response.isin(["yes", "no"])
        ]
        cohort_label = f"{ex_condition} + {ex_treatment} + {ex_sample}"
    else:
        cohort = data[
            (data.condition == "melanoma")
            & (data.treatment == "miraclib")
            & (data.sample_type == "PBMC")
            & data.response.isin(["yes", "no"])
        ]
        cohort_label = "melanoma + miraclib + PBMC (required cohort)"

    st.markdown(f"**Cohort:** {cohort_label}")

    if cohort.empty or cohort.response.nunique() < 2:
        st.warning("Not enough data in this cohort to compare responders vs non-responders.")
    else:
        c1, c2, c3 = st.columns(3)
        c1.metric("Samples", cohort.sample_id.nunique())
        c2.metric("Responders", cohort[cohort.response == "yes"].subject_id.nunique())
        c3.metric("Non-responders", cohort[cohort.response == "no"].subject_id.nunique())

        results = run_response_tests(cohort)

        st.subheader("Boxplots: % of total by response × timepoint")
        populations = sorted(cohort.population.unique())
        timepoints_in_cohort = sorted(cohort.timepoint.unique())

        # 5 populations as rows, timepoints as columns of one big plotly figure
        # via facet_col on timepoint and facet_row on population
        fig = px.box(
            cohort,
            x="response", y="percentage",
            color="response",
            facet_row="population", facet_col="timepoint",
            points="outliers",
            category_orders={
                "response": ["no", "yes"],
                "population": populations,
                "timepoint": timepoints_in_cohort,
            },
            color_discrete_map={"no": "#d9534f", "yes": "#3a7ca5"},
            height=180 * len(populations),
        )
        fig.update_yaxes(matches=None, title="% of total")
        fig.update_xaxes(title=None)
        # Strip the "population=" / "timepoint=" prefixes from facet labels.
        fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
        fig.update_layout(margin=dict(t=40, b=10), showlegend=False)
        st.plotly_chart(fig, use_container_width=True)

        st.subheader("Statistical results")
        st.caption("Mann-Whitney U (two-sided), Benjamini-Hochberg FDR across all tests.")
        display = results.copy()
        # Bold-highlight significant rows using a style function
        def highlight_sig(row):
            return [
                "background-color: #d4edda; color: #155724" if row.significant_q05 else ""
                for _ in row
            ]
        st.dataframe(
            display.style.apply(highlight_sig, axis=1).format({
                "median_yes": "{:.3f}",
                "median_no": "{:.3f}",
                "p_value": "{:.4f}",
                "q_value_bh": "{:.4f}",
            }),
            use_container_width=True,
        )

        # Interpretation block: pull live numbers from the results table
        n_sig = int(results.significant_q05.sum())
        lowest = results.sort_values("q_value_bh").iloc[0]
        st.subheader("Interpretation")
        if n_sig == 0:
            st.info(
                f"**Primary finding:** No immune cell population showed a "
                f"statistically significant difference between responders and "
                f"non-responders after FDR correction across "
                f"{len(results)} comparisons. Lowest q-value: "
                f"**{lowest.q_value_bh:.3f}** "
                f"({lowest.population} @ t={int(lowest.timepoint)}).\n\n"
                "**Conclusion:** In this cohort, bulk cell-type proportions do "
                "not appear to be a useful predictive biomarker for response. "
                "Finer immune phenotyping or other modalities would be needed "
                "to identify response predictors."
            )
            nominal = results[(results.p_value < 0.05) & ~results.significant_q05]
            if len(nominal):
                st.warning(
                    f"**Exploratory observations (hypothesis-generating only):** "
                    f"{len(nominal)} comparison(s) had uncorrected p < 0.05 but "
                    "did not survive multiple-testing correction. These are listed "
                    "for completeness and should not be interpreted as predictive "
                    "without replication in an independent cohort."
                )
                st.dataframe(
                    nominal[["population", "timepoint", "median_yes", "median_no", "p_value", "q_value_bh"]],
                    use_container_width=True,
                )
        else:
            sig = results[results.significant_q05]
            st.success(
                f"**Primary finding:** {n_sig} of {len(results)} comparisons "
                "showed a statistically significant difference (FDR q < 0.05)."
            )
            st.dataframe(sig, use_container_width=True)


# Tab 4: Baseline subset (Part 4)
with tab_subset:
    st.header("Baseline subset analysis")
    st.caption(
        "Part 4 deliverable. Baseline cohort definition is fixed by the spec: "
        "**melanoma + miraclib + PBMC, time_from_treatment_start = 0**."
    )

    baseline = data[
        (data.condition == "melanoma")
        & (data.treatment == "miraclib")
        & (data.sample_type == "PBMC")
        & (data.timepoint == 0)
    ]

    c1, c2, c3 = st.columns(3)
    c1.metric("Samples", baseline.sample_id.nunique())
    c2.metric("Subjects", baseline.subject_id.nunique())
    c3.metric("Projects represented", baseline.project_id.nunique())

    st.divider()

    col_a, col_b, col_c = st.columns(3)

    with col_a:
        st.subheader("Samples per project")
        proj = (
            baseline.drop_duplicates("sample_id")
            .project_id.value_counts().reset_index()
        )
        proj.columns = ["project_id", "samples"]
        fig = px.bar(proj.sort_values("project_id"), x="project_id", y="samples", text="samples")
        fig.update_layout(showlegend=False, height=300, margin=dict(t=10, b=10))
        st.plotly_chart(fig, use_container_width=True)

    with col_b:
        st.subheader("Subjects by response")
        sub = baseline.drop_duplicates("subject_id")
        resp = sub.response.value_counts().reset_index()
        resp.columns = ["response", "subjects"]
        fig = px.bar(
            resp, x="response", y="subjects", color="response", text="subjects",
            color_discrete_map={"yes": "#3a7ca5", "no": "#d9534f"},
        )
        fig.update_layout(showlegend=False, height=300, margin=dict(t=10, b=10))
        st.plotly_chart(fig, use_container_width=True)

    with col_c:
        st.subheader("Subjects by sex")
        sex_c = sub.sex.value_counts().reset_index()
        sex_c.columns = ["sex", "subjects"]
        fig = px.bar(
            sex_c, x="sex", y="subjects", color="sex", text="subjects",
            color_discrete_map={"M": "#3a7ca5", "F": "#d96a8a"},
        )
        fig.update_layout(showlegend=False, height=300, margin=dict(t=10, b=10))
        st.plotly_chart(fig, use_container_width=True)

    st.divider()
    st.subheader("Key metric: B-cell count, melanoma males who responded at t=0")

    male_resp_bcells = baseline[
        (baseline.sex == "M")
        & (baseline.response == "yes")
        & (baseline.population == "b_cell")
    ]
    n = len(male_resp_bcells)
    mean_b = male_resp_bcells["count"].mean()
    median_b = male_resp_bcells["count"].median()
    std_b = male_resp_bcells["count"].std()

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("n samples", n)
    c2.metric("Mean B-cells", f"{mean_b:.2f}")
    c3.metric("Median B-cells", f"{median_b:.0f}")
    c4.metric("Std dev", f"{std_b:.0f}")

    fig = px.histogram(
        male_resp_bcells, x="count", nbins=25,
        title=None, color_discrete_sequence=["#3a7ca5"],
    )
    fig.add_vline(x=mean_b, line_dash="dash", line_color="black",
                  annotation_text=f"mean = {mean_b:.2f}", annotation_position="top")
    fig.update_layout(height=320, margin=dict(t=10, b=10),
                      xaxis_title="B-cell count", yaxis_title="samples")
    st.plotly_chart(fig, use_container_width=True)
