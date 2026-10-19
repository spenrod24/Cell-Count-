""" Initialize and populate the clinical-trial SQLite database. Run from the repository root: python load_data.py """
import csv
import sqlite3
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CSV_PATH = ROOT / 'cell-count.csv'
DB_PATH = ROOT / 'clinical_trial.db'
CELL_TYPES = ('b_cell', 'cd8_t_cell', 'cd4_t_cell', 'nk_cell', 'monocyte')

SCHEMA = '''
PRAGMA foreign_keys = ON;
CREATE TABLE projects (project_id TEXT PRIMARY KEY);
CREATE TABLE subjects (
    project_id TEXT NOT NULL REFERENCES projects(project_id),
    subject_id TEXT NOT NULL,
    condition TEXT NOT NULL,
    age INTEGER NOT NULL CHECK(age >= 0),
    sex TEXT NOT NULL,
    treatment TEXT NOT NULL,
    response TEXT NOT NULL,
    PRIMARY KEY (project_id, subject_id)
);
CREATE TABLE samples (
    sample_id TEXT PRIMARY KEY,
    project_id TEXT NOT NULL,
    subject_id TEXT NOT NULL,
    sample_type TEXT NOT NULL,
    time_from_treatment_start INTEGER NOT NULL,
    FOREIGN KEY(project_id, subject_id) REFERENCES subjects(project_id, subject_id)
);
CREATE TABLE cell_counts (
    sample_id TEXT NOT NULL REFERENCES samples(sample_id),
    cell_type TEXT NOT NULL CHECK(cell_type IN ('b_cell','cd8_t_cell','cd4_t_cell','nk_cell','monocyte')),
    cell_count INTEGER NOT NULL CHECK(cell_count >= 0),
    PRIMARY KEY(sample_id, cell_type)
);
CREATE INDEX idx_samples_subject_time ON samples(project_id, subject_id, time_from_treatment_start);
CREATE INDEX idx_subjects_treatment_response ON subjects(treatment, response);
CREATE INDEX idx_cell_counts_type ON cell_counts(cell_type);
'''

def load():
    if not CSV_PATH.is_file():
        raise FileNotFoundError(f'Missing input CSV: {CSV_PATH}')

    # Read and validate before replacing any existing database.
    with CSV_PATH.open(newline='', encoding='utf-8-sig') as source:
        reader = csv.DictReader(source)
        expected = {'project','subject','condition','age','sex','treatment','response',
                    'sample','sample_type','time_from_treatment_start',*CELL_TYPES}
        if not expected.issubset(set(reader.fieldnames or [])):
            raise ValueError(f'CSV missing columns: {sorted(expected - set(reader.fieldnames or []))}')
        rows = list(reader)
    temp_path = ROOT / 'clinical_trial.loading.db'
    if temp_path.exists():
        temp_path.unlink()
    try:
        with sqlite3.connect(temp_path) as conn:
            conn.executescript(SCHEMA)
            for row in rows:
                project = row['project'].strip()
                subject = row['subject'].strip()
                conn.execute('INSERT OR IGNORE INTO projects VALUES (?)', (project,))
                demographics = (project, subject, row['condition'], int(row['age']),
                                row['sex'], row['treatment'], row['response'])
                conn.execute('INSERT OR IGNORE INTO subjects VALUES (?,?,?,?,?,?,?)', demographics)
                saved = conn.execute('SELECT project_id,subject_id,condition,age,sex,treatment,response FROM subjects WHERE project_id=? AND subject_id=?', (project, subject)).fetchone()
                if saved != demographics:
                    raise ValueError(f'Conflicting subject attributes: {project}/{subject}')
                conn.execute('INSERT INTO samples VALUES (?,?,?,?,?)',
                             (row['sample'], project, subject, row['sample_type'], int(row['time_from_treatment_start'])))
                conn.executemany('INSERT INTO cell_counts VALUES (?,?,?)',
                                 [(row['sample'], cell_type, int(row[cell_type])) for cell_type in CELL_TYPES])
        temp_path.replace(DB_PATH)
    finally:
        if temp_path.exists():
            temp_path.unlink()
    print(f'Loaded {len(rows):,} samples and {len(rows)*len(CELL_TYPES):,} cell counts into {DB_PATH.name}')

if __name__ == '__main__':
    load()
