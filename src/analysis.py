"""
Part 2: Frequency summary table

For every (sample, population) pair, we compute:
  total_count: the sample's total across all 5 populations
  count: the count for that population
  percentage: count / total_count * 100

Output: outputs/frequency_table.csv with columns
  sample, total_count, population, count, percentage
"""

from __future__ import annotations

import sqlite3
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DB_PATH = ROOT / "loblaw.db"
OUT_PATH = ROOT / "outputs" / "frequency_table.csv"


# One query does whole computation. Window function gives each sample's
# total, and we keep the long-format shape
FREQUENCY_QUERY = """
SELECT
    sample_id   AS sample,
    SUM(count) OVER (PARTITION BY sample_id) AS total_count,
    population,
    count,
    ROUND(100.0 * count / SUM(count) OVER (PARTITION BY sample_id), 4) AS percentage
FROM cell_counts
ORDER BY sample_id, population;
"""

def compute_frequencies(conn: sqlite3.Connection) -> pd.DataFrame:
    """Return the Part 2 frequency table as a DataFrame."""
    return pd.read_sql_query(FREQUENCY_QUERY, conn)

def main() -> None:
    if not DB_PATH.exists():
        raise SystemExit(
            f"ERROR: {DB_PATH} not found. Run `python load_data.py` first."
        )

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

    print("Computing frequency table ...")
    with sqlite3.connect(DB_PATH) as conn:
        df = compute_frequencies(conn)

    # Check: every sample's 5 percentages sums to 100
    sums = df.groupby("sample")["percentage"].sum()
    max_dev = (sums - 100).abs().max()
    if max_dev > 0.01:
        raise SystemExit(
            f"ERROR: percentages don't sum to 100 for some samples "
            f"(max deviation: {max_dev:.4f})"
        )

    df.to_csv(OUT_PATH, index=False)
    print(f"  {len(df):,} rows written to {OUT_PATH.relative_to(ROOT)}")
    print(f"  Samples × populations: {df['sample'].nunique()} × {df['population'].nunique()}")
    print(f"  Per-sample percentage sum check: max |sum - 100| = {max_dev:.2e}")

    print("\nFirst 10 rows:")
    print(df.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
