# Loblaw Bio — Cell Count Analysis

**Live dashboard:** https://loblaw-bio-analysis-pessccwsk7xugbg4bgmhcn.streamlit.app/

## How to run

In a GitHub Codespace (or any terminal with `make`, `python3`, and `pip`):

```bash
make setup       # install dependencies from requirements.txt
make pipeline    # run all four parts end-to-end
make dashboard   # launch the Streamlit dashboard locally
```

`make pipeline` initializes `loblaw.db` from `cell-count.csv` and writes every output (frequency table, statistical results, boxplots, subset CSVs) into `outputs/`. Total runtime: ~10 seconds.

## Database schema

Four tables in 3NF:

**`projects`**

| Column | Type | Notes |
|---|---|---|
| `project_id` | TEXT | Primary key |

**`subjects`** — one row per study subject

| Column | Type | Notes |
|---|---|---|
| `subject_id` | TEXT | Primary key |
| `project_id` | TEXT | Foreign key → `projects` |
| `condition` | TEXT | melanoma / carcinoma / healthy |
| `age` | INTEGER | |
| `sex` | TEXT | M / F |

**`samples`** — one row per biological sample (a subject typically has 3, one per timepoint)

| Column | Type | Notes |
|---|---|---|
| `sample_id` | TEXT | Primary key |
| `subject_id` | TEXT | Foreign key → `subjects` |
| `sample_type` | TEXT | PBMC / WB |
| `treatment` | TEXT | miraclib / phauximab / none |
| `response` | TEXT | yes / no / NULL (for untreated) |
| `time_from_treatment_start` | INTEGER | days |

**`cell_counts`** — one row per (sample, population) pair, long format

| Column | Type | Notes |
|---|---|---|
| `id` | INTEGER | Primary key |
| `sample_id` | TEXT | Foreign key → `samples` |
| `population` | TEXT | b_cell / cd8_t_cell / cd4_t_cell / nk_cell / monocyte |
| `count` | INTEGER | |

UNIQUE constraint on `(sample_id, population)`. Relationships: `projects` 1→N `subjects` 1→N `samples` 1→5 `cell_counts`.

Rationale:
Subject demographics concern subjects, not samples. The unprocessed CSV duplicates the condition, age, and sex for every one of a subject's three timepoint rows. This denormalized structure wastes space. Dividing a table into smaller tables ensures that each fact can only be stated once.
The cell_counts table is in a long format (one row per sample × population). This is the most important scaling decision: introducing a new sample (regulatory T cell, NKT cell, or a new activation marker) would now be a simple INSERT, instead of an ALTER TABLE along with a migration for every affected query.
Even though it is only an ID for now, a separate projects table was created because, in production, a project can include things like a PI, IRB, a start date, and project status. Adding these later will be an expensive exercise.
All foreign keys have individual indexes and composite index on the exact filter combinations that Parts 3 and 4 use (treatment, time_from_treatment_start, and sample_type).

Scaling: At 100 projects, 1000 samples, 5 populations = 500K cell_counts rows. SQLite handles this without problem, with sub-100ms queries for specific data. The same structure in Postgres works with only the edits of INTEGER PRIMARY KEY AUTOINCREMENT to SERIAL. For a more demanding workload with hundreds of thousands (or millions) of samples, real-time dashboards, and multiple users, placing a star schema in a columnar warehouse (like DuckDB, BigQuery, or Snowflake) on top of the same OLTP base would improve query latency without impacting ingestion. The long-format structure of the cell_counts table is the main contributor. A wide format would lead to schema migrations every time an assay panel was modified.

## Code structure

```
.
├── load_data.py              # Part 1: schema + CSV → SQLite (root, executable)
├── src/
│   ├── analysis.py           # Part 2: per-sample frequency table
│   ├── statistics.py         # Part 3: responder vs non-responder stats
│   └── subset_analysis.py    # Part 4: baseline subset queries
├── dashboard/app.py          # Streamlit dashboard
├── Makefile                  # setup / pipeline / dashboard targets
└── requirements.txt
```

One script per assignment part. Each script can be independently run and tested like how a bioinformatics group would segment ingestion, summarization, statistics, and ad-hoc analysis.

A couple notes on design choices:
The frequency calculation (Part 2) is written in SQL, not pandas. A window function (SUM(count) OVER (PARTITION BY sample_id)) is a single query to compute per-sample totals — the database handles heavy lifting, not just data storage.
Part 4 implements a CTE in order to write the baseline-cohort filter once and reference it in the rest of the four sub-queries. If the cohort definition changes, this is the only edit required.
Part 3 implements Mann-Whitney U with Benjamini-Hochberg FDR correction to the five population subsets at three time points for a total of fifteen tests, with rank-biserial effect sizes and associated p-values. This approach is supported by an EDA notebook (notebooks/eda_stats_check.ipynb) where Shapiro-Wilk test of normality is not passed in 19 of the 30 groups, and where Levene's test of equality of variances is not passed in 5 of the 15 test comparisons, thus violating the assumptions of t-tests.
The dashboard computes the statistics from Part 3 in real-time, as opposed to relying on an already generated CSV, so the "explore other cohorts" option runs the same analysis on the subsets selected by the user.

## Part 3 results

After FDR correction applied to all 15 comparisons, there were no statistically significant differences in the relative frequency of any of the immune cell populations between those who responded to miraclib and those who did not respond to miraclib, based on melanoma PBMC samples (lowest q-value: 0.22; |effect size r| ≤ 0.11). The PBMC cell-type proportions thus do not appear to be a reasonable candidate predictive biomarker for the response to miraclib for this cohort.

## Part 4 results

Baseline cohort (melanoma + miraclib + PBMC, t=0): 656 samples from 656 subjects (prj1: 384, prj3: 272, prj2: 0). 331 responders / 325 non-responders; 344 male / 312 female. Melanoma males at t=0 had a mean B-cell count of 10401.28.
