# Loblaw Bio Immune Cell Analysis

This repo contains a Python analysis pipeline and dashboard for immune cell population data from a clinical trial:

1. Load `cell-count.csv` into a SQLite database.
2. Calculate relative immune cell population frequencies per sample.
3. Compare melanoma PBMC miraclib responders and non-responders.
4. Generate summary outputs and plots.
5. Provide an interactive dashboard.

## Running the project

Install dependencies:

```bash
make setup

```
