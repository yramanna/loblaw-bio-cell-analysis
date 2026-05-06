import sqlite3
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent
CSV_PATH = ROOT / "cell-count.csv"
DB_PATH = ROOT / "cell_counts.db"

CELL_POPULATIONS = [
    "b_cell",
    "cd8_t_cell",
    "cd4_t_cell",
    "nk_cell",
    "monocyte",
]


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Normalize the actual CSV columns to names used internally.

    Actual CSV columns include:
    project, subject, condition, age, sex, treatment, response, sample,
    sample_type, time_from_treatment_start, and the five cell populations.
    """
    df = df.copy()
    df.columns = [c.strip() for c in df.columns]

    rename_map = {
        "sample": "sample_id",
        "subject": "patient_id",
        "condition": "indication",
        "sex": "gender",
    }

    df = df.rename(columns=rename_map)

    required = {
        "project",
        "patient_id",
        "indication",
        "age",
        "gender",
        "treatment",
        "response",
        "sample_id",
        "sample_type",
        "time_from_treatment_start",
        *CELL_POPULATIONS,
    }

    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in cell-count.csv: {missing}")

    return df


def initialize_database(conn: sqlite3.Connection) -> None:
    cur = conn.cursor()

    cur.executescript(
        """
        DROP TABLE IF EXISTS cell_counts;
        DROP TABLE IF EXISTS samples;
        DROP TABLE IF EXISTS patients;
        DROP TABLE IF EXISTS projects;

        CREATE TABLE projects (
            project_id TEXT PRIMARY KEY
        );

        CREATE TABLE patients (
            patient_id TEXT PRIMARY KEY,
            indication TEXT NOT NULL,
            age INTEGER,
            treatment TEXT NOT NULL,
            response TEXT NOT NULL,
            gender TEXT NOT NULL
        );

        CREATE TABLE samples (
            sample_id TEXT PRIMARY KEY,
            patient_id TEXT NOT NULL,
            project_id TEXT NOT NULL,
            sample_type TEXT NOT NULL,
            time_from_treatment_start REAL NOT NULL,
            FOREIGN KEY (patient_id) REFERENCES patients(patient_id),
            FOREIGN KEY (project_id) REFERENCES projects(project_id)
        );

        CREATE TABLE cell_counts (
            sample_id TEXT NOT NULL,
            population TEXT NOT NULL,
            count INTEGER NOT NULL,
            PRIMARY KEY (sample_id, population),
            FOREIGN KEY (sample_id) REFERENCES samples(sample_id)
        );
        """
    )

    conn.commit()


def load_data(conn: sqlite3.Connection, df: pd.DataFrame) -> None:
    cur = conn.cursor()

    for _, row in df.iterrows():
        project_id = str(row["project"])
        patient_id = str(row["patient_id"])
        sample_id = str(row["sample_id"])

        cur.execute(
            """
            INSERT OR IGNORE INTO projects (project_id)
            VALUES (?)
            """,
            (project_id,),
        )

        cur.execute(
            """
            INSERT OR IGNORE INTO patients
                (patient_id, indication, age, treatment, response, gender)
            VALUES (?, ?, ?, ?, ?, ?)
            """,
            (
                patient_id,
                str(row["indication"]),
                int(row["age"]),
                str(row["treatment"]),
                str(row["response"]),
                str(row["gender"]),
            ),
        )

        cur.execute(
            """
            INSERT OR REPLACE INTO samples
                (
                    sample_id,
                    patient_id,
                    project_id,
                    sample_type,
                    time_from_treatment_start
                )
            VALUES (?, ?, ?, ?, ?)
            """,
            (
                sample_id,
                patient_id,
                project_id,
                str(row["sample_type"]),
                float(row["time_from_treatment_start"]),
            ),
        )

        for population in CELL_POPULATIONS:
            cur.execute(
                """
                INSERT OR REPLACE INTO cell_counts
                    (sample_id, population, count)
                VALUES (?, ?, ?)
                """,
                (
                    sample_id,
                    population,
                    int(row[population]),
                ),
            )

    conn.commit()


def main() -> None:
    if not CSV_PATH.exists():
        raise FileNotFoundError(
            "cell-count.csv was not found in the repository root."
        )

    df = pd.read_csv(CSV_PATH)
    df = normalize_columns(df)

    conn = sqlite3.connect(DB_PATH)
    try:
        initialize_database(conn)
        load_data(conn, df)
    finally:
        conn.close()

    print(f"Database created: {DB_PATH}")


if __name__ == "__main__":
    main()
