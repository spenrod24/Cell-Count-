# Loblaw Bio — Clinical Trial Immune Cell Analysis

SQLite-backed Python pipeline and interactive Streamlit dashboard for the provided `cell-count.csv`.

## Run in GitHub Codespaces

1. Open this repository in **Code → Codespaces → Create codespace on main**.
2. In the terminal, run:

   ```bash
   make setup
   make pipeline
   make dashboard
   ```

3. When Codespaces detects port **8501**, open the forwarded port in your browser. If it does not appear automatically, open the **Ports** panel and select **Open in Browser** for port 8501. Keep `make dashboard` running.

**Dashboard link:** [Open the locally running dashboard](http://localhost:8501). In Codespaces, use the forwarded port URL displayed in the Ports panel; a permanent hosted URL has not been configured.

`make setup` installs dependencies; `make pipeline` recreates `clinical_trial.db` from the CSV and regenerates all result files; `make dashboard` starts Streamlit. The targets must be invoked from the repository root. Python 3.10+ is recommended (GitHub Codespaces Python 3.11+ works).

## Project structure

- `cell-count.csv` — original input data
- `load_data.py` — creates normalized SQLite tables `projects`, `subjects`, `samples`, `cell_counts`
- `analyze_cells.py` — outputs `cell_frequencies.csv` (five rows per sample)
- `analyze_response.py` — outputs sample-level and subject-level frequency CSVs, `response_statistics.csv`, and `response_boxplots.png`
- `analyze_baseline.py` — outputs `baseline_samples.csv`, `baseline_by_project.csv`, `baseline_by_response.csv`, and `baseline_by_sex.csv`
- `dashboard.py` — interactive SQLite-backed tabs for Parts 2–4
- `requirements.txt`, `Makefile` — reproducible environment and execution

## Analysis methods

**Part 2:** Each population's frequency is `100 * cell_count / sum(all five counts for the sample)`. Zero-total samples have `NULL` frequency in the database-backed dashboard; the exported Part 2 CSV uses `0.0` for these samples.

**Part 3:** Filter melanoma, miraclib, PBMC, and response yes/no. Average repeated sample frequencies for each unique `(project, subject)` before comparing response groups using two-sided Mann–Whitney U tests. Adjust five p-values with Benjamini–Hochberg false-discovery-rate correction; q < 0.05 is the prespecified exploratory threshold. The boxplots and dashboard reflect subject-level averages. These are associations, not validated predictive biomarkers.

**Part 4:** Filter melanoma, miraclib, PBMC, time_from_treatment_start = 0. Count samples per project and distinct subjects by response and sex.

## Troubleshooting

- `No module named streamlit`: run `make setup`.
- Missing SQLite database: run `make pipeline` before `make dashboard`.
- Port 8501 already in use: stop the prior dashboard process before restarting.
