import sqlite3
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import mannwhitneyu


ROOT = Path(__file__).resolve().parent
DB_PATH = ROOT / "cell_counts.db"
OUTPUT_DIR = ROOT / "outputs"
OUTPUT_DIR.mkdir(exist_ok=True)


def connect() -> sqlite3.Connection:
    if not DB_PATH.exists():
        raise FileNotFoundError(
            "Database not found. Run `python load_data.py` first."
        )
    return sqlite3.connect(DB_PATH)


def create_summary_table(conn: sqlite3.Connection) -> pd.DataFrame:
    query = """
    WITH totals AS (
        SELECT
            sample_id,
            SUM(count) AS total_count
        FROM cell_counts
        GROUP BY sample_id
    )
    SELECT
        cc.sample_id AS sample,
        totals.total_count,
        cc.population,
        cc.count,
        ROUND(100.0 * cc.count / totals.total_count, 4) AS percentage
    FROM cell_counts cc
    JOIN totals
        ON cc.sample_id = totals.sample_id
    ORDER BY cc.sample_id, cc.population;
    """
    return pd.read_sql_query(query, conn)


def miraclib_response_data(conn: sqlite3.Connection) -> pd.DataFrame:
    query = """
    WITH totals AS (
        SELECT
            sample_id,
            SUM(count) AS total_count
        FROM cell_counts
        GROUP BY sample_id
    )
    SELECT
        s.sample_id AS sample,
        p.patient_id,
        p.indication,
        p.age,
        p.treatment,
        p.response,
        p.gender,
        s.project_id AS project,
        s.sample_type,
        s.time_from_treatment_start,
        cc.population,
        cc.count,
        100.0 * cc.count / totals.total_count AS percentage
    FROM samples s
    JOIN patients p
        ON s.patient_id = p.patient_id
    JOIN cell_counts cc
        ON s.sample_id = cc.sample_id
    JOIN totals
        ON s.sample_id = totals.sample_id
    WHERE LOWER(p.indication) = 'melanoma'
      AND LOWER(p.treatment) = 'miraclib'
      AND LOWER(s.sample_type) = 'pbmc'
      AND LOWER(p.response) IN ('yes', 'no');
    """
    return pd.read_sql_query(query, conn)


def run_statistics(df: pd.DataFrame) -> pd.DataFrame:
    results = []

    for population in sorted(df["population"].unique()):
        subset = df[df["population"] == population]

        responders = subset[subset["response"].str.lower() == "yes"]["percentage"]
        non_responders = subset[subset["response"].str.lower() == "no"]["percentage"]

        if len(responders) == 0 or len(non_responders) == 0:
            p_value = None
            statistic = None
            significant = False
        else:
            test = mannwhitneyu(
                responders,
                non_responders,
                alternative="two-sided",
            )
            statistic = test.statistic
            p_value = test.pvalue
            significant = p_value < 0.05

        results.append(
            {
                "population": population,
                "n_responder_samples": len(responders),
                "n_non_responder_samples": len(non_responders),
                "mean_percent_responders": responders.mean(),
                "mean_percent_non_responders": non_responders.mean(),
                "mannwhitneyu_statistic": statistic,
                "p_value": p_value,
                "significant_at_0.05": significant,
            }
        )

    return pd.DataFrame(results)


def make_boxplot(df: pd.DataFrame) -> None:
    populations = sorted(df["population"].unique())
    responses = ["yes", "no"]

    data = []
    labels = []

    for population in populations:
        for response in responses:
            values = df[
                (df["population"] == population)
                & (df["response"].str.lower() == response)
            ]["percentage"]
            data.append(values)
            labels.append(f"{population}\n{response}")

    plt.figure(figsize=(12, 6))
    plt.boxplot(data, labels=labels)
    plt.ylabel("Relative frequency (%)")
    plt.xlabel("Population and response group")
    plt.title("Melanoma PBMC miraclib samples: responders vs non-responders")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    out_path = OUTPUT_DIR / "miraclib_response_boxplot.png"
    plt.savefig(out_path, dpi=200)
    plt.close()


def subset_analysis(conn: sqlite3.Connection) -> dict:
    baseline_query = """
    SELECT
        s.sample_id AS sample,
        s.project_id AS project,
        p.patient_id,
        p.indication,
        p.age,
        p.treatment,
        p.response,
        p.gender,
        s.sample_type,
        s.time_from_treatment_start
    FROM samples s
    JOIN patients p
        ON s.patient_id = p.patient_id
    WHERE LOWER(p.indication) = 'melanoma'
      AND LOWER(p.treatment) = 'miraclib'
      AND LOWER(s.sample_type) = 'pbmc'
      AND s.time_from_treatment_start = 0;
    """

    baseline = pd.read_sql_query(baseline_query, conn)

    project_counts = (
        baseline.groupby("project")["sample"]
        .nunique()
        .reset_index(name="sample_count")
    )

    response_counts = (
        baseline.groupby("response")["patient_id"]
        .nunique()
        .reset_index(name="subject_count")
    )

    gender_counts = (
        baseline.groupby("gender")["patient_id"]
        .nunique()
        .reset_index(name="subject_count")
    )

    avg_b_cell_query = """
    SELECT
        AVG(cc.count) AS avg_b_cells
    FROM samples s
    JOIN patients p
        ON s.patient_id = p.patient_id
    JOIN cell_counts cc
        ON s.sample_id = cc.sample_id
    WHERE LOWER(p.indication) = 'melanoma'
      AND LOWER(p.treatment) = 'miraclib'
      AND LOWER(s.sample_type) = 'pbmc'
      AND s.time_from_treatment_start = 0
      AND LOWER(p.gender) = 'm'
      AND LOWER(p.response) = 'yes'
      AND cc.population = 'b_cell';
    """

    avg_b_cells = pd.read_sql_query(avg_b_cell_query, conn)
    avg_value = avg_b_cells.loc[0, "avg_b_cells"]

    return {
        "baseline_samples": baseline,
        "project_counts": project_counts,
        "response_counts": response_counts,
        "gender_counts": gender_counts,
        "avg_b_cells_melanoma_male_responders_time0": avg_value,
    }


def main() -> None:
    with connect() as conn:
        summary = create_summary_table(conn)
        summary.to_csv(OUTPUT_DIR / "summary_relative_frequencies.csv", index=False)

        response_df = miraclib_response_data(conn)
        response_df.to_csv(OUTPUT_DIR / "miraclib_response_data.csv", index=False)

        stats = run_statistics(response_df)
        stats.to_csv(OUTPUT_DIR / "miraclib_response_statistics.csv", index=False)

        if not response_df.empty:
            make_boxplot(response_df)

        subsets = subset_analysis(conn)
        subsets["baseline_samples"].to_csv(
            OUTPUT_DIR / "baseline_melanoma_pbmc_miraclib_samples.csv",
            index=False,
        )
        subsets["project_counts"].to_csv(
            OUTPUT_DIR / "baseline_project_counts.csv",
            index=False,
        )
        subsets["response_counts"].to_csv(
            OUTPUT_DIR / "baseline_response_counts.csv",
            index=False,
        )
        subsets["gender_counts"].to_csv(
            OUTPUT_DIR / "baseline_gender_counts.csv",
            index=False,
        )

        avg_b = subsets["avg_b_cells_melanoma_male_responders_time0"]
        with open(OUTPUT_DIR / "avg_b_cells_male_responders_time0.txt", "w") as f:
            if pd.isna(avg_b):
                f.write("No matching samples found.\n")
            else:
                f.write(f"{avg_b:.2f}\n")

    print("Pipeline complete. Outputs written to outputs/")


if __name__ == "__main__":
    main()
