"""Fig. 4 SGA – step 1: parse raw colony data, reproduce the original top-200 list,
compute rdh54-specific genetic-interaction scores and GO-category tests.
Run:  python code/11_sga_scores.py  (from any directory)
Inputs : raw_data/SGA_screen/data_heatmap_completo_replica123_sergio_log2fold.csv  (original, unmodified)
         reference/go-basic.obo.gz, reference/pombase.gaf.gz
Outputs: output/sga_long_72h.csv, output/top200_reproduced.csv,
         output/sga_gene_scores.csv, output/sga_enrichment.csv
"""
import os as _os; _os.chdir(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))  # run from the Reproducibility folder
import re, gzip, os
import numpy as np, pandas as pd, networkx as nx
from scipy import stats
os.makedirs('output/SGA_screen', exist_ok=True)
RAW = 'raw_data/SGA_screen/data_heatmap_completo_replica123_sergio_log2fold.csv'
d = pd.read_csv(RAW, sep=';')
# label example: 10_72hYESHU4601_rep2-_B_Batch1  -> plate 10, 72h, YES+HU, query 4601, replicate B
rx = re.compile(r'^(\d+)_(\d+h)(YESHU|YES)(\d+)(?:_rep\d)?-_([ABC])_')
p = d.label.str.extract(rx); p.columns = ['plate', 'time', 'medium', 'query', 'rep']
d = pd.concat([d, p], axis=1); assert d['query'].notna().all()
d['v'] = d.log2fold_value_normalized2          # NA kept as NA (= colony absent)
# query 4601 = rdh54Δ ; query 3206 = control query ; YESHU = YES + 15 mM hydroxyurea

# ---------- 1. Reproduce the ORIGINAL ranking (Top_200_SGA.R): NA->0, 72h YESHU 4601 only
o = d.copy(); o['x'] = o.log2fold_value_normalized2.fillna(0)
o = o[(o.time == '72h') & (o.medium == 'YESHU') & (o['query'] == '4601')]
tops = {r: set(o[o.rep == r].sort_values('x').head(500).strain) for r in 'ABC'}
keep = (tops['A'] & tops['B']) | (tops['A'] & tops['C']) | (tops['B'] & tops['C'])
t = o[o.strain.isin(keep)].groupby('strain').x.sum() / 3
t = t.sort_values().head(200)
names = d.drop_duplicates('strain').set_index('strain').gene2.astype(str)
pd.DataFrame({'gene': names.reindex(t.index), 'score_original': t, 'rank': range(1, len(t) + 1)}).to_csv('output/SGA_screen/top200_reproduced.csv')

# ---------- 2. Interaction scores (this re-analysis)
d72 = d[d.time == '72h']; d72.to_csv('output/SGA_screen/sga_long_72h.csv', index=False)
rep = d72.groupby(['strain', 'query', 'medium', 'rep']).v.mean().unstack(['query', 'medium', 'rep'])
G = pd.DataFrame(index=rep.index); G['gene'] = names.reindex(G.index)
for q in ('4601', '3206'):
    for m in ('YES', 'YESHU'):
        sub = rep.xs((q, m), axis=1, level=(0, 1))
        G[f'{q}_{m}_mean'] = sub.mean(axis=1); G[f'{q}_{m}_n'] = sub.notna().sum(axis=1)
G['score_HU'] = G['4601_YESHU_mean'] - G['3206_YESHU_mean']
G['score_YES'] = G['4601_YES_mean'] - G['3206_YES_mean']
G['evaluable_HU'] = (G['4601_YESHU_n'] >= 2) & (G['3206_YESHU_n'] >= 1)
G['absent_in_rdh54_HU'] = (G['4601_YESHU_n'] == 0) & (G['3206_YESHU_n'] >= 1)

# ---------- 3. GO categories (PomBase GAF, propagated through is_a / part_of)
gr = nx.DiGraph(); cur = None
for line in gzip.open('raw_data/SGA_screen/annotation/go-basic.obo.gz', 'rt'):
    line = line.strip()
    if line == '[Term]': cur = None
    elif line.startswith('id: GO:'): cur = line[4:]
    elif cur and line.startswith('is_a: '): gr.add_edge(cur, line[6:16])
    elif cur and line.startswith('relationship: part_of '): gr.add_edge(cur, line[22:32])
ann = {}
with gzip.open('raw_data/SGA_screen/annotation/pombase.gaf.gz', 'rt') as f:
    for l in f:
        if l.startswith('!'): continue
        x = l.split('\t')
        if 'NOT' in x[3]: continue
        ann.setdefault(x[1], set()).add(x[4])
anc = {}
def ancestors(t):
    if t not in anc: anc[t] = ({t} | nx.descendants(gr, t)) if t in gr else {t}
    return anc[t]
full = {g: set().union(*[ancestors(t) for t in ts]) for g, ts in ann.items()}
SETS = {'DNA repair': ['GO:0006281'], 'DNA replication': ['GO:0006260'],
        'Actin cytoskeleton': ['GO:0030036', 'GO:0015629', 'GO:0003779', 'GO:0007015'],
        'Microtubule cytoskeleton': ['GO:0000226', 'GO:0015630', 'GO:0008017'],
        'Myosin complex': ['GO:0016459']}
for s, ts in SETS.items(): G[s] = [bool(full.get(g, set()) & set(ts)) for g in G.index]
E = G[G.evaluable_HU].copy(); N = len(E); G['rank_HU'] = E.score_HU.rank(method='first')
top = E.nsmallest(200, 'score_HU'); res = []
for s in SETS:
    ins = E[s]; k = int(top[s].sum()); K = int(ins.sum())
    res.append(dict(category=s, genes_in_screen=K, in_top200=k, expected=round(K * 200 / N, 1),
        fold=round(k / (K * 200 / N), 2) if K else np.nan,
        P_hypergeom_top200=stats.hypergeom.sf(k - 1, N, K, 200),
        P_rank_HU=stats.mannwhitneyu(E.loc[ins, 'score_HU'], E.loc[~ins, 'score_HU'], alternative='less').pvalue,
        P_rank_YES=stats.mannwhitneyu(E.loc[ins, 'score_YES'].dropna(), E.loc[~ins, 'score_YES'].dropna(), alternative='less').pvalue,
        median_score_HU=E.loc[ins, 'score_HU'].median()))
G.to_csv('output/SGA_screen/sga_gene_scores.csv'); pd.DataFrame(res).to_csv('output/SGA_screen/sga_enrichment.csv', index=False)
print('evaluable deletions (HU):', N); print(pd.DataFrame(res).round(4).to_string(index=False))
