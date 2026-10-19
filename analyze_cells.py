"""Calculate per-sample immune cell relative frequencies.

Run: python analyze_cells.py
Creates cell_frequencies.csv and prints the first 15 results.
Requires clinical_trial.db created by python load_data.py.
"""
import csv
import sqlite3
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DB_PATH = ROOT / 'clinical_trial.db'
OUTPUT_PATH = ROOT / 'cell_frequencies.csv'
COLUMNS = ('sample', 'total_count', 'population', 'count', 'percentage')

QUERY = '''
WITH counts AS (
    SELECT sample_id, cell_type, cell_count,
           SUM(cell_count) OVER (PARTITION BY sample_id) AS total_count
    FROM cell_counts
)
SELECT sample_id AS sample, total_count, cell_type AS population,
       cell_count AS count,
       CASE WHEN total_count = 0 THEN 0.0
            ELSE ROUND(100.0 * cell_count / total_count, 6) END AS percentage
FROM counts
ORDER BY sample_id,
    CASE cell_type
        WHEN 'b_cell' THEN 1
        WHEN 'cd8_t_cell' THEN 2
        WHEN 'cd4_t_cell' THEN 3
        WHEN 'nk_cell' THEN 4
        WHEN 'monocyte' THEN 5
    END;
'''

def get_frequencies(connection):
    """Return one row per sample/population with percentages in [0, 100]."""
    return connection.execute(QUERY)

def main():
    if not DB_PATH.is_file():
        raise FileNotFoundError('Missing clinical_trial.db; run python load_data.py first.')
    with sqlite3.connect(DB_PATH) as connection:
        cursor = get_frequencies(connection)
        with OUTPUT_PATH.open('w', newline='', encoding='utf-8') as output:
            writer = csv.writer(output)
            writer.writerow(COLUMNS)
            preview = []
            row_count = 0
            for row in cursor:
                writer.writerow(row)
                if len(preview) < 15:
                    preview.append(row)
                row_count += 1
    print(', '.join(COLUMNS))
    for row in preview:
        print(', '.join(map(str, row)))
    print(f'\nSaved {row_count:,} rows to {OUTPUT_PATH.name}')

if __name__ == '__main__':
    main()
