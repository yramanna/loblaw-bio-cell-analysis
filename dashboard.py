from pathlib import Path
import sqlite3

import pandas as pd
import plotly.express as px
import streamlit as st
from scipy.stats import mannwhitneyu


ROOT = Path(__file__).resolve().parent
DB_PATH = ROOT / "cell_counts.db"


@st.cache_data
def load_summary() -> pd.DataFrame:
    conn = sqlite3.connect(DB_PATH)

    query = """
    WITH totals AS (
        SELECT sample_id, SUM(count) AS total_count
        FROM cell_counts
        GROUP BY sample_id
    )
    SELECT
        cc.sample_id AS sample,
        totals.total_count,
        cc.population,
        cc.count,
        100.0 * cc.count / totals.total_count AS percentage,
        p.patient_id,
        p.indication,
        p.age,
        p.treatment,
        p.response,
        p.gender,
        s.project_id AS project,
        s.sample_type,
        s.time_from_treatment_start
    FROM cell_counts cc
    JOIN totals
        ON cc.sample_id = totals.sample_id
    JOIN samples s
        ON cc.sample_id = s.sample_id
    JOIN patients p
        ON s.patient_id = p.patient_id;
    """

    df = pd.read_sql_query(query, conn)
    conn.close()
    return df


def compute_stats(df: pd.DataFrame) -> pd.DataFrame:
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
            test = mannwhitneyu(responders, non_responders, alternative="two-sided")
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


def main() -> None:
    st.set_page_config(
        page_title="Loblaw Bio Immune Cell Analysis",
        layout="wide",
    )

    st.title("Loblaw Bio Immune Cell Population Analysis")

    if not DB_PATH.exists():
        st.error("Database not found. Run `make pipeline` first.")
        return

    df = load_summary()

    st.sidebar.header("Filters")

    indications = st.sidebar.multiselect(
        "Condition / indication",
        sorted(df["indication"].dropna().unique()),
        default=sorted(df["indication"].dropna().unique()),
    )

    treatments = st.sidebar.multiselect(
        "Treatment",
        sorted(df["treatment"].dropna().unique()),
        default=sorted(df["treatment"].dropna().unique()),
    )

    sample_types = st.sidebar.multiselect(
        "Sample type",
        sorted(df["sample_type"].dropna().unique()),
        default=sorted(df["sample_type"].dropna().unique()),
    )

    populations = st.sidebar.multiselect(
        "Cell population",
        sorted(df["population"].dropna().unique()),
        default=sorted(df["population"].dropna().unique()),
    )

    filtered = df[
        df["indication"].isin(indications)
        & df["treatment"].isin(treatments)
        & df["sample_type"].isin(sample_types)
        & df["population"].isin(populations)
    ]

    st.header("Part 2: Relative Frequency Summary")
    st.dataframe(
        filtered[
            [
                "sample",
                "total_count",
                "population",
                "count",
                "percentage",
                "project",
                "patient_id",
                "indication",
                "treatment",
                "response",
                "gender",
                "sample_type",
                "time_from_treatment_start",
            ]
        ],
        use_container_width=True,
    )

    st.header("Part 3: Miraclib Responder vs Non-Responder Analysis")

    response_df = df[
        (df["indication"].str.lower() == "melanoma")
        & (df["treatment"].str.lower() == "miraclib")
        & (df["sample_type"].str.lower() == "pbmc")
        & (df["response"].str.lower().isin(["yes", "no"]))
    ]

    if response_df.empty:
        st.warning("No melanoma PBMC miraclib responder/non-responder samples found.")
    else:
        fig = px.box(
            response_df,
            x="population",
            y="percentage",
            color="response",
            points="all",
            title="Relative frequency by response status",
            labels={
                "population": "Cell population",
                "percentage": "Relative frequency (%)",
                "response": "Response",
            },
        )
        st.plotly_chart(fig, use_container_width=True)

        stats = compute_stats(response_df)
        st.subheader("Statistical test results")
        st.write(
            "Mann-Whitney U tests compare relative frequencies between responders and non-responders."
        )
        st.dataframe(stats, use_container_width=True)

        significant = stats[stats["significant_at_0.05"]]
        if significant.empty:
            st.info("No cell populations were significant at p < 0.05.")
        else:
            st.success(
                "Significant populations at p < 0.05: "
                + ", ".join(significant["population"].tolist())
            )

    st.header("Part 4: Baseline Melanoma PBMC Miraclib Samples")

    baseline = df[
        (df["indication"].str.lower() == "melanoma")
        & (df["treatment"].str.lower() == "miraclib")
        & (df["sample_type"].str.lower() == "pbmc")
        & (df["time_from_treatment_start"] == 0)
    ]

    baseline_samples = baseline.drop_duplicates("sample")

    st.subheader("Matching baseline samples")
    st.dataframe(
        baseline_samples[
            [
                "sample",
                "project",
                "patient_id",
                "indication",
                "treatment",
                "response",
                "gender",
                "sample_type",
                "time_from_treatment_start",
            ]
        ],
        use_container_width=True,
    )

    col1, col2, col3 = st.columns(3)

    with col1:
        st.subheader("Samples per project")
        st.dataframe(
            baseline_samples.groupby("project")["sample"]
            .nunique()
            .reset_index(name="sample_count"),
            use_container_width=True,
        )

    with col2:
        st.subheader("Subjects by response")
        st.dataframe(
            baseline_samples.groupby("response")["patient_id"]
            .nunique()
            .reset_index(name="subject_count"),
            use_container_width=True,
        )

    with col3:
        st.subheader("Subjects by gender")
        st.dataframe(
            baseline_samples.groupby("gender")["patient_id"]
            .nunique()
            .reset_index(name="subject_count"),
            use_container_width=True,
        )

    b_cell_subset = baseline[
        (baseline["population"] == "b_cell")
        & (baseline["response"].str.lower() == "yes")
        & (baseline["gender"].str.lower() == "m")
    ]

    avg_b_cells = b_cell_subset["count"].mean()

    st.subheader("Average B-cell count")
    if pd.isna(avg_b_cells):
        st.write("No matching melanoma male responder samples at time 0 were found.")
    else:
        st.metric(
            "Melanoma male responders, miraclib, PBMC, time 0",
            f"{avg_b_cells:.2f}",
        )

    st.caption(
        "The current dataset includes miraclib, phauximab, and none. "
        "The schema can support future treatment arms or AI model comparisons, including quintazide."
    )


if __name__ == "__main__":
    main()
