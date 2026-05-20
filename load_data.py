"""
Part 1: Data Management

Reads cell-count.csv, creates a normalized SQLite schema, and loads all rows

Schema (4 tables, 3NF):
  projects      one row per project
  subjects      one row per subject  (FK -> projects)
  samples       one row per sample   (FK -> subjects)
  cell_counts   one row per (sample, population) pair  (FK -> samples)

Long-format cell_counts table
"""

from __future__ import annotations

import sqlite3
import sys
from pathlib import Path

import pandas as pd

# Paths resolved relative to this file
ROOT = Path(__file__).resolve().parent
CSV_PATH = ROOT / "cell-count.csv"
DB_PATH = ROOT / "loblaw.db"

# The five immune cell populations measured in this study
POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]


SCHEMA_SQL = """
PRAGMA foreign_keys = ON;

DROP TABLE IF EXISTS cell_counts;
DROP TABLE IF EXISTS samples;
DROP TABLE IF EXISTS subjects;
DROP TABLE IF EXISTS projects;

CREATE TABLE projects (
    project_id  TEXT PRIMARY KEY,
    created_at  TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE subjects (
    subject_id  TEXT PRIMARY KEY,
    project_id  TEXT NOT NULL,
    condition   TEXT NOT NULL,
    age         INTEGER,
    sex         TEXT CHECK (sex IN ('M', 'F')),
    FOREIGN KEY (project_id) REFERENCES projects(project_id)
);

CREATE TABLE samples (
    sample_id                   TEXT PRIMARY KEY,
    subject_id                  TEXT NOT NULL,
    sample_type                 TEXT NOT NULL,
    treatment                   TEXT,
    response                    TEXT CHECK (response IN ('yes', 'no') OR response IS NULL),
    time_from_treatment_start   INTEGER,
    FOREIGN KEY (subject_id) REFERENCES subjects(subject_id)
);

CREATE TABLE cell_counts (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    sample_id   TEXT NOT NULL,
    population  TEXT NOT NULL,
    count       INTEGER NOT NULL CHECK (count >= 0),
    FOREIGN KEY (sample_id) REFERENCES samples(sample_id),
    UNIQUE (sample_id, population)
);

-- Indexes for the filters used in Parts 3 & 4.
CREATE INDEX idx_subjects_project    ON subjects(project_id);
CREATE INDEX idx_subjects_condition  ON subjects(condition);
CREATE INDEX idx_samples_subject     ON samples(subject_id);
CREATE INDEX idx_samples_treatment   ON samples(treatment);
CREATE INDEX idx_samples_time        ON samples(time_from_treatment_start);
CREATE INDEX idx_samples_filter      ON samples(treatment, time_from_treatment_start, sample_type);
CREATE INDEX idx_cellcounts_sample   ON cell_counts(sample_id);
CREATE INDEX idx_cellcounts_pop      ON cell_counts(population);
"""


def load() -> None:
    """Build the database from scratch and load all rows from cell-count.csv."""

    if not CSV_PATH.exists():
        sys.exit(f"ERROR: {CSV_PATH} not found. Place cell-count.csv in the repo root.")

    print(f"Reading {CSV_PATH.name} ...")
    df = pd.read_csv(CSV_PATH)
    print(f"  {len(df):,} rows, {df.shape[1]} columns")

    # Validate expected columns are present before touching DB
    expected = {
        "project", "subject", "condition", "age", "sex",
        "treatment", "response", "sample", "sample_type",
        "time_from_treatment_start", *POPULATIONS,
    }
    missing = expected - set(df.columns)
    if missing:
        sys.exit(f"ERROR: CSV is missing required columns: {sorted(missing)}")

    # Normalize 'response': empty/NaN -> None so SQLite stores NULL
    df["response"] = df["response"].where(df["response"].notna(), None)

    # Remove any existing DB file so script is idempotent
    if DB_PATH.exists():
        DB_PATH.unlink()
        print(f"  Removed existing {DB_PATH.name}")

    print(f"Creating database {DB_PATH.name} ...")
    with sqlite3.connect(DB_PATH) as conn:
        conn.executescript(SCHEMA_SQL)

        # projects: one row per unique project_id
        projects = df[["project"]].drop_duplicates().rename(columns={"project": "project_id"})
        projects.to_sql("projects", conn, if_exists="append", index=False)
        print(f"  projects:    {len(projects):>7,} rows")

        # subjects: one row per subject. Demographics are constant per subject
        # drop_duplicates collapses the 3 timepoint rows into 1 subject row
        subjects = (
            df[["subject", "project", "condition", "age", "sex"]]
            .drop_duplicates(subset=["subject"])
            .rename(columns={"subject": "subject_id", "project": "project_id"})
        )
        subjects.to_sql("subjects", conn, if_exists="append", index=False)
        print(f"  subjects:    {len(subjects):>7,} rows")

        # samples: one row per sample_id
        samples = df[[
            "sample", "subject", "sample_type", "treatment",
            "response", "time_from_treatment_start",
        ]].rename(columns={"sample": "sample_id", "subject": "subject_id"})
        samples.to_sql("samples", conn, if_exists="append", index=False)
        print(f"  samples:     {len(samples):>7,} rows")

        # cell_counts: long format. Melt the 5 population columns into rows
        cell_counts = df.melt(
            id_vars=["sample"],
            value_vars=POPULATIONS,
            var_name="population",
            value_name="count",
        ).rename(columns={"sample": "sample_id"})
        cell_counts.to_sql("cell_counts", conn, if_exists="append", index=False)
        print(f"  cell_counts: {len(cell_counts):>7,} rows")

        # Check: cell_counts should be exactly samples * populations.
        expected_counts = len(samples) * len(POPULATIONS)
        if len(cell_counts) != expected_counts:
            sys.exit(
                f"ERROR: expected {expected_counts:,} cell_counts rows, "
                f"got {len(cell_counts):,}"
            )

        conn.commit()

    print(f"\nDone. Database written to {DB_PATH}")


if __name__ == "__main__":
    load()
