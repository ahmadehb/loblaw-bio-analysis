"""
Part 4: Baseline subset analysis

Defines the baseline cohort:
  melanoma + miraclib + PBMC + time_from_treatment_start = 0

Reports:
  1. Samples per project
  2. Responders vs non-responders (subject-level counts)
  3. Males vs females (subject-level counts)
  4. Average B-cell count for melanoma males who responded at t=0  

Outputs CSVs in outputs/ and prints a summary block to stdout
"""

from __future__ import annotations

import sqlite3
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DB_PATH = ROOT / "loblaw.db"
OUT_DIR = ROOT / "outputs"


# Reusable CTE defines the baseline cohort, used by all queries below
# so filter is written exactly once
BASELINE_CTE = """
WITH baseline AS (
    SELECT
        sub.subject_id, sub.project_id, sub.sex, sub.condition,
        s.sample_id, s.treatment, s.response, s.time_from_treatment_start
    FROM subjects sub
    JOIN samples s ON s.subject_id = sub.subject_id
    WHERE sub.condition = 'melanoma'
      AND s.treatment   = 'miraclib'
      AND s.sample_type = 'PBMC'
      AND s.time_from_treatment_start = 0
)
"""


def samples_per_project(conn: sqlite3.Connection) -> pd.DataFrame:
    """Q4.1: how many baseline samples per project."""
    return pd.read_sql_query(
        BASELINE_CTE + """
        SELECT project_id, COUNT(DISTINCT sample_id) AS sample_count
        FROM baseline
        GROUP BY project_id
        ORDER BY project_id;
        """,
        conn,
    )


def subjects_by_response(conn: sqlite3.Connection) -> pd.DataFrame:
    """Q4.2: responder vs non-responder subject counts."""
    return pd.read_sql_query(
        BASELINE_CTE + """
        SELECT response, COUNT(DISTINCT subject_id) AS subject_count
        FROM baseline
        GROUP BY response
        ORDER BY response;
        """,
        conn,
    )


def subjects_by_sex(conn: sqlite3.Connection) -> pd.DataFrame:
    """Q4.3: male vs female subject counts."""
    return pd.read_sql_query(
        BASELINE_CTE + """
        SELECT sex, COUNT(DISTINCT subject_id) AS subject_count
        FROM baseline
        GROUP BY sex
        ORDER BY sex;
        """,
        conn,
    )


def mean_bcells_male_responders(conn: sqlite3.Connection) -> float:
    """
    Q4.4: average B-cell count for melanoma males who responded at t=0.

    Uses raw cell counts (not percentages) — the spec asks for "average number
    of B cells".
    """
    row = pd.read_sql_query(
        BASELINE_CTE + """
        SELECT AVG(cc.count) AS mean_bcells, COUNT(*) AS n_samples
        FROM baseline b
        JOIN cell_counts cc ON cc.sample_id = b.sample_id
        WHERE b.sex = 'M'
          AND b.response = 'yes'
          AND cc.population = 'b_cell';
        """,
        conn,
    ).iloc[0]
    return float(row["mean_bcells"]), int(row["n_samples"])


def main() -> None:
    if not DB_PATH.exists():
        raise SystemExit(f"ERROR: {DB_PATH} not found. Run `python load_data.py` first.")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Querying baseline cohort: melanoma + miraclib + PBMC + t=0 ...")
    with sqlite3.connect(DB_PATH) as conn:
        # Cohort size for the header.
        total = pd.read_sql_query(
            BASELINE_CTE + "SELECT COUNT(DISTINCT sample_id) AS n FROM baseline;",
            conn,
        ).iloc[0]["n"]
        print(f"  Baseline cohort: {total} samples\n")

        proj = samples_per_project(conn)
        resp = subjects_by_response(conn)
        sex  = subjects_by_sex(conn)
        mean_b, n_b = mean_bcells_male_responders(conn)

    # Persist each breakdown table.
    proj.to_csv(OUT_DIR / "subset_samples_per_project.csv", index=False)
    resp.to_csv(OUT_DIR / "subset_subjects_by_response.csv", index=False)
    sex.to_csv(OUT_DIR / "subset_subjects_by_sex.csv", index=False)

    male_responder_summary = pd.DataFrame([{
        "cohort": "melanoma_male_responders_t0",
        "population": "b_cell",
        "n_samples": n_b,
        "mean_count": round(mean_b, 2),
    }])
    male_responder_summary.to_csv(
        OUT_DIR / "subset_mean_bcells_male_responders.csv", index=False,
    )

    # Pretty print to stdout.
    print("Samples per project:")
    print(proj.to_string(index=False))
    print("\nSubjects by response (responder=yes, non-responder=no):")
    print(resp.to_string(index=False))
    print("\nSubjects by sex:")
    print(sex.to_string(index=False))
    print(f"\nAverage B-cell count, melanoma males who responded at t=0:")
    print(f"  n = {n_b} samples")
    print(f"  mean = {mean_b:.2f}")

    print(f"\nWrote 4 CSVs to {OUT_DIR.relative_to(ROOT)}/")


if __name__ == "__main__":
    main()
