# Loblaw Bio Immune Cell Analysis

This project analyzes immune cell counts from `cell-count.csv`, stores the data in SQLite, generates the required analysis outputs, and provides a Streamlit dashboard.

## Run in GitHub Codespaces

Install dependencies:

    make setup

Run the full pipeline:

    make pipeline

This creates `cell_counts.db` and writes output tables/plots to `outputs/`.

Start the dashboard:

    make dashboard

Dashboard link:

    http://localhost:8501

## Database schema

The SQLite database has four tables:

- `projects(project_id)`
- `patients(patient_id, indication, age, treatment, response, gender)`
- `samples(sample_id, patient_id, project_id, sample_type, time_from_treatment_start)`
- `cell_counts(sample_id, population, count)`

The input CSV uses `subject`, `condition`, `sex`, and `sample`. These are loaded as `patient_id`, `indication`, `gender`, and `sample_id`.

This schema separates project, patient, sample, and cell-count information. It avoids repeating metadata for every cell population and supports multiple samples per subject. It also scales well to hundreds of projects and thousands of samples because new samples and cell populations can be added as rows without changing the table structure.

## Code structure

- `load_data.py`: creates the SQLite database and loads `cell-count.csv`
- `analysis.py`: calculates relative frequencies, runs statistical analysis, and generates output files
- `dashboard.py`: Streamlit dashboard for exploring results
- `Makefile`: provides `setup`, `pipeline`, and `dashboard` commands
- `requirements.txt`: Python dependencies

The code is split this way so loading, analysis, and visualization can be run independently while still using the same SQLite database.
