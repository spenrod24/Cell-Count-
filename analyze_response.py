"""Part 3: exploratory PBMC response comparison, using subject-level averages.

Run: python analyze_response.py
Requires: python load_data.py (and clinical_trial.db).
Dependencies: numpy, pandas, scipy, matplotlib.
"""
from pathlib import Path
import sqlite3
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
DB = ROOT / 'clinical_trial.db'
CELL_TYPES = ['b_cell', 'cd8_t_cell', 'cd4_t_cell', 'nk_cell', 'monocyte']

SQL = '''
WITH freq AS (
 SELECT s.sample_id AS sample, s.project_id, s.subject_id,
        s.time_from_treatment_start AS timepoint,
        sub.response, c.cell_type AS population,
        c.cell_count AS count,
        SUM(c.cell_count) OVER (PARTITION BY s.sample_id) AS total_count
 FROM samples s
 JOIN subjects sub ON sub.project_id=s.project_id AND sub.subject_id=s.subject_id
 JOIN cell_counts c ON c.sample_id=s.sample_id
 WHERE LOWER(TRIM(sub.condition))='melanoma'
   AND LOWER(TRIM(sub.treatment))='miraclib'
   AND UPPER(TRIM(s.sample_type))='PBMC'
   AND LOWER(TRIM(sub.response)) IN ('yes','no')
)
SELECT *, CASE WHEN total_count=0 THEN NULL
 ELSE 100.0*count/total_count END AS percentage FROM freq
'''

def bh_adjust(pvalues):
    p = np.asarray(pvalues, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.minimum.accumulate((ranked * len(p) / np.arange(1, len(p)+1))[::-1])[::-1]
    out = np.empty(len(p), dtype=float)
    out[order] = np.clip(adjusted, 0, 1)
    return out

def run():
    if not DB.exists():
        raise FileNotFoundError('Run python load_data.py first.')
    with sqlite3.connect(DB) as conn:
        samples = pd.read_sql_query(SQL, conn)
    if samples.empty:
        raise ValueError('No eligible melanoma/miraclib/PBMC samples.')
    samples.to_csv(ROOT / 'response_sample_frequencies.csv', index=False)
    # Each subject contributes one average per cell type, avoiding repeated-measure pseudoreplication.
    subjects = (samples.dropna(subset=['percentage'])
                .groupby(['project_id','subject_id','response','population'], as_index=False)
                .agg(percentage=('percentage','mean'), n_samples=('sample','nunique')))
    subjects.to_csv(ROOT / 'response_subject_frequencies.csv', index=False)
    records = []
    for population in CELL_TYPES:
        subset = subjects[subjects.population == population]
        yes = subset[subset.response == 'yes'].percentage.to_numpy()
        no = subset[subset.response == 'no'].percentage.to_numpy()
        if not len(yes) or not len(no):
            raise ValueError(f'Insufficient response groups for {population}')
        test = mannwhitneyu(yes, no, alternative='two-sided', method='asymptotic')
        records.append({'population':population,'n_responders':len(yes),
                        'n_nonresponders':len(no),
                        'responder_median_pct':float(np.median(yes)),
                        'nonresponder_median_pct':float(np.median(no)),
                        'median_difference_pp':float(np.median(yes)-np.median(no)),
                        'p_value':float(test.pvalue)})
    results = pd.DataFrame(records)
    results['q_value_bh'] = bh_adjust(results.p_value)
    results['significant_fdr_0_05'] = results.q_value_bh < 0.05
    results.to_csv(ROOT / 'response_statistics.csv', index=False, float_format='%.9g')

    fig, axes = plt.subplots(1, 5, figsize=(17, 5), sharey=True)
    for ax, pop in zip(axes, CELL_TYPES):
        part = subjects[subjects.population == pop]
        a = part[part.response=='yes'].percentage.to_numpy()
        b = part[part.response=='no'].percentage.to_numpy()
        ax.boxplot([a,b], tick_labels=[f'Yes\n(n={len(a)})', f'No\n(n={len(b)})'],
                   showfliers=False, patch_artist=True,
                   boxprops={'facecolor':'#b5d5e8'}, medianprops={'color':'#152f46'})
        ax.set_title(pop.replace('_',' ')+'\nq='+format(results.loc[results.population==pop,'q_value_bh'].iloc[0],'.3g'))
        ax.grid(axis='y', alpha=.25)
    axes[0].set_ylabel('Cell frequency (%) — subject mean')
    fig.suptitle('Melanoma · miraclib · PBMC | Responders vs non-responders', fontsize=13)
    fig.text(.5,.01,'Subject-level mean across available timepoints; two-sided Mann–Whitney U; BH-adjusted q values',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.04,1,.91))
    fig.savefig(ROOT/'response_boxplots.png', dpi=180)
    plt.close(fig)
    print(f"Eligible samples: {samples['sample'].nunique():,}; subjects: {subjects[['project_id','subject_id']].drop_duplicates().shape[0]:,}")
    print(results.to_string(index=False, float_format=lambda x:f'{x:.5g}'))
    print('Saved response_sample_frequencies.csv, response_subject_frequencies.csv, response_statistics.csv, response_boxplots.png')
    return results

if __name__ == '__main__':
    run()
