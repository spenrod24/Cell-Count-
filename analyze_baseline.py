"""Part 4: melanoma / miraclib / PBMC baseline subset and summaries.

Run: python analyze_baseline.py (after python load_data.py).
"""
import csv
import sqlite3
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DB = ROOT / 'clinical_trial.db'

BASELINE_SQL = '''
SELECT s.sample_id AS sample, s.project_id AS project,
       s.subject_id AS subject, sub.condition, sub.treatment,
       sub.response, sub.sex, s.sample_type,
       s.time_from_treatment_start
FROM samples AS s
JOIN subjects AS sub
  ON sub.project_id = s.project_id AND sub.subject_id = s.subject_id
WHERE LOWER(TRIM(sub.condition)) = 'melanoma'
  AND LOWER(TRIM(sub.treatment)) = 'miraclib'
  AND UPPER(TRIM(s.sample_type)) = 'PBMC'
  AND s.time_from_treatment_start = 0
ORDER BY s.project_id, s.subject_id, s.sample_id
'''

PROJECT_SQL = '''
SELECT project, COUNT(*) AS sample_count, COUNT(DISTINCT subject) AS subject_count
FROM baseline GROUP BY project ORDER BY project
'''

# Group at subject granularity to avoid counting multiple baseline samples twice.
SUBJECTS_SQL = '''
SELECT project, subject, response, sex FROM baseline
GROUP BY project, subject, response, sex
'''


def export_csv(path, rows, columns):
    with path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def print_table(title, rows, columns):
    print('\n' + title)
    print(' | '.join(columns))
    for row in rows:
        print(' | '.join(str(row[col]) for col in columns))


def run():
    if not DB.is_file():
        raise FileNotFoundError('Database not found. Run python load_data.py first.')
    with sqlite3.connect(DB) as conn:
        conn.row_factory = sqlite3.Row
        baseline = [dict(row) for row in conn.execute(BASELINE_SQL)]
        conn.execute('CREATE TEMP VIEW baseline AS ' + BASELINE_SQL.replace('ORDER BY s.project_id, s.subject_id, s.sample_id', ''))
        projects = [dict(row) for row in conn.execute(PROJECT_SQL)]
        response = [dict(row) for row in conn.execute('''
            SELECT LOWER(TRIM(response)) AS response, COUNT(*) AS subject_count
            FROM (''' + SUBJECTS_SQL + ''') GROUP BY LOWER(TRIM(response)) ORDER BY response
        ''')]
        sex = [dict(row) for row in conn.execute('''
            SELECT UPPER(TRIM(sex)) AS sex, COUNT(*) AS subject_count
            FROM (''' + SUBJECTS_SQL + ''') GROUP BY UPPER(TRIM(sex)) ORDER BY sex
        ''')]
    export_csv(ROOT / 'baseline_samples.csv', baseline,
               ['sample','project','subject','condition','treatment','response','sex','sample_type','time_from_treatment_start'])
    export_csv(ROOT / 'baseline_by_project.csv', projects, ['project','sample_count','subject_count'])
    export_csv(ROOT / 'baseline_by_response.csv', response, ['response','subject_count'])
    export_csv(ROOT / 'baseline_by_sex.csv', sex, ['sex','subject_count'])
    print(f'Eligible baseline samples: {len(baseline)}')
    print_table('Samples by project', projects, ['project','sample_count','subject_count'])
    print_table('Subjects by response', response, ['response','subject_count'])
    print_table('Subjects by sex', sex, ['sex','subject_count'])
    print('\nSaved baseline_samples.csv, baseline_by_project.csv, baseline_by_response.csv, baseline_by_sex.csv')
    return baseline, projects, response, sex


if __name__ == '__main__':
    run()
