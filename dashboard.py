"""Interactive SQLite-backed clinical trial dashboard. Run: streamlit run dashboard.py"""
from pathlib import Path
import sqlite3
import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st
from scipy.stats import mannwhitneyu
from analyze_response import bh_adjust, CELL_TYPES

ROOT = Path(__file__).resolve().parent
DB = ROOT / 'clinical_trial.db'
st.set_page_config(page_title='Loblaw Bio | Immune Cell Explorer', layout='wide')
st.title('Loblaw Bio · Immune Cell Explorer')
st.caption('Exploratory clinical-trial analytics | Frequencies, response comparison, and baseline cohorts')
if not DB.exists():
    st.error('Database missing. Run `make pipeline` first.')
    st.stop()

@st.cache_data

def load_data(db_mtime):
    sql = '''WITH raw AS (
      SELECT s.sample_id AS sample, s.project_id AS project, s.subject_id AS subject,
      sub.condition, sub.treatment, sub.response, sub.sex, sub.age,
      s.sample_type, s.time_from_treatment_start AS timepoint,
      c.cell_type AS population, c.cell_count AS count,
      SUM(c.cell_count) OVER (PARTITION BY s.sample_id) AS total_count
      FROM samples s JOIN subjects sub ON sub.project_id=s.project_id AND sub.subject_id=s.subject_id
      JOIN cell_counts c ON c.sample_id=s.sample_id)
      SELECT *, CASE WHEN total_count=0 THEN NULL ELSE 100.0*count/total_count END AS percentage FROM raw'''
    with sqlite3.connect(DB) as conn:
        return pd.read_sql_query(sql, conn)

df = load_data(DB.stat().st_mtime_ns)
frequency_tab, response_tab, baseline_tab = st.tabs(['Cell frequencies · Part 2', 'Response analysis · Part 3', 'Baseline subset · Part 4'])
with frequency_tab:
    st.subheader('Relative frequency by sample')
    a,b,c=st.columns(3)
    project=a.multiselect('Project',sorted(df.project.unique()),default=sorted(df.project.unique()),key='freq_proj')
    condition=b.multiselect('Condition',sorted(df.condition.unique()),default=sorted(df.condition.unique()))
    sample_type=c.multiselect('Sample type',sorted(df.sample_type.unique()),default=sorted(df.sample_type.unique()))
    subset=df[df.project.isin(project)&df.condition.isin(condition)&df.sample_type.isin(sample_type)]
    sample_ids=sorted(subset['sample'].unique())
    chosen=st.selectbox('Inspect a sample',sample_ids) if sample_ids else None
    if chosen:
        one=subset[subset['sample']==chosen]
        st.metric('Total cells in selected sample',f'{int(one.total_count.iloc[0]):,}')
        st.bar_chart(one.set_index('population')['percentage'],y_label='Frequency (%)')
    st.caption('Complete filtered table; download as CSV for all matching samples.')
    table=subset[['sample','total_count','population','count','percentage']].copy()
    st.dataframe(table,hide_index=True,use_container_width=True,height=330)
    st.download_button('Download filtered frequencies',table.to_csv(index=False).encode(),'cell_frequencies_filtered.csv','text/csv')
with response_tab:
    st.subheader('Melanoma · miraclib · PBMC: responders vs non-responders')
    st.caption('Patient-level means across timepoints avoid treating repeated samples as independent. Exploratory, not a validated prediction model.')
    eligible=df[(df.condition.str.lower()=='melanoma')&(df.treatment.str.lower()=='miraclib')&(df.sample_type.str.upper()=='PBMC')&(df.response.str.lower().isin(['yes','no']))]
    projects=st.multiselect('Projects to compare',sorted(eligible.project.unique()),default=sorted(eligible.project.unique()),key='response_proj')
    eligible=eligible[eligible.project.isin(projects)]
    subject=(eligible.dropna(subset=['percentage']).groupby(['project','subject','response','population'],as_index=False)['percentage'].mean())
    if subject.empty:
        st.info('No matching subjects.')
    else:
        st.metric('Subjects',len(subject[['project','subject']].drop_duplicates()))
        fig,axes=plt.subplots(1,5,figsize=(16,4.5),sharey=True)
        records=[]
        for ax,pop in zip(axes,CELL_TYPES):
            group=subject[subject.population==pop]
            yes=group[group.response.str.lower()=='yes'].percentage.to_numpy()
            no=group[group.response.str.lower()=='no'].percentage.to_numpy()
            if len(yes) and len(no):
                ax.boxplot([yes,no],tick_labels=['Yes','No'],showfliers=False)
                p=float(mannwhitneyu(yes,no,alternative='two-sided',method='asymptotic').pvalue)
                records.append({'population':pop,'responders':len(yes),'non_responders':len(no),'responder_median_pct':float(pd.Series(yes).median()),'non_responder_median_pct':float(pd.Series(no).median()),'p_value':p})
            ax.set_title(pop.replace('_',' '))
            ax.grid(axis='y',alpha=.2)
        axes[0].set_ylabel('Frequency (%) · patient mean')
        fig.tight_layout()
        st.pyplot(fig)
        plt.close(fig)
        if records:
            stats=pd.DataFrame(records)
            stats['q_value_bh']=bh_adjust(stats.p_value)
            stats['significant_q_lt_0_05']=stats.q_value_bh<0.05
            st.dataframe(stats,hide_index=True,use_container_width=True)
            st.download_button('Download comparison statistics',stats.to_csv(index=False).encode(),'response_statistics_filtered.csv','text/csv')
            st.info('Two-sided Mann–Whitney U tests; Benjamini–Hochberg adjustment across tested populations. Associations do not establish predictive performance.')
with baseline_tab:
    st.subheader('Baseline: melanoma · miraclib · PBMC · time = 0')
    base=df[(df.condition.str.lower()=='melanoma')&(df.treatment.str.lower()=='miraclib')&(df.sample_type.str.upper()=='PBMC')&(df.timepoint==0)]
    base=base.drop_duplicates(subset=['sample'])
    project_filter=st.multiselect('Project filter',sorted(base.project.unique()),default=sorted(base.project.unique()),key='baseline_proj')
    base=base[base.project.isin(project_filter)]
    st.metric('Baseline samples',len(base))
    subjects=base.drop_duplicates(subset=['project','subject'])
    left,middle,right=st.columns(3)
    with left:
        st.markdown('**Samples by project**')
        st.dataframe(base.groupby('project').agg(samples=('sample','nunique'),subjects=('subject','nunique')).reset_index(),hide_index=True)
    with middle:
        st.markdown('**Unique subjects by response**')
        st.dataframe(subjects.groupby('response').size().rename('subjects').reset_index(),hide_index=True)
    with right:
        st.markdown('**Unique subjects by sex**')
        st.dataframe(subjects.groupby('sex').size().rename('subjects').reset_index(),hide_index=True)
    st.markdown('**Qualifying samples**')
    st.dataframe(base[['sample','project','subject','condition','treatment','response','sex','sample_type','timepoint']],hide_index=True,use_container_width=True)
    st.download_button('Download baseline subset',base.to_csv(index=False).encode(),'baseline_filtered.csv','text/csv')
