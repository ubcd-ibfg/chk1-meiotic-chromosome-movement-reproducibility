"""Fig. 4 SGA – step 2: panels A–C and the rank plot.  Run after 11_sga_scores.py.
python code/12_sga_panels.py  -> output/SGA_screen/panels/, output/SGA_screen/source_data_*.csv
"""
import os as _os; _os.chdir(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))  # run from the Reproducibility folder
import os, logging, numpy as np, pandas as pd, matplotlib
logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
plt.rcParams.update({'font.family': ['Arial', 'Liberation Sans'], 'font.size': 7, 'axes.spines.top': False,
    'axes.spines.right': False, 'pdf.fonttype': 42, 'svg.fonttype': 'none', 'axes.linewidth': 0.7,
    'xtick.major.width': 0.7, 'ytick.major.width': 0.7})
O = 'output/SGA_screen/panels/'; os.makedirs(O, exist_ok=True)
G = pd.read_csv('output/SGA_screen/sga_gene_scores.csv', index_col=0)
d = pd.read_csv('output/SGA_screen/sga_long_72h.csv', dtype={'query': str})
ACT, MYO, REP, SPB, GREY, HU = '#D6338A', '#1B9E77', '#2B6CB0', '#7B61C4', '#CFCFCF', '#C8901A'
# genes highlighted (category from PomBase GO, see 11_sga_scores.py)
HL = {  # gene: (colour, label offset in points)
    'rho4': (ACT, (8, -2)), 'rdi1': (ACT, (8, -4)), 'ssm4': (ACT, (8, 3)), 'cap1': (ACT, (8, 0)),
    'mid1': (ACT, (8, 2)), 'myo51': (MYO, (-34, 4)), 'msh2': (REP, (-30, -2)), 'hus1': (REP, (-26, -8)),
    'mcp6': (SPB, (5, -9)), 'bqt1': (SPB, (-26, -4))}
def it(g): return f'${g}\\Delta$'
def save(f, n):
    for e in ('pdf', 'svg'): f.savefig(O + n + '.' + e, bbox_inches='tight', dpi=600)
    f.savefig(O + n + '.png', dpi=300, bbox_inches='tight'); plt.close(f)

# ---------------- A  scheme
f, ax = plt.subplots(figsize=(7.0, 1.55)); ax.axis('off'); ax.set_xlim(-1.5, 102); ax.set_ylim(0, 22)
def box(x, y, w, h, t, c='#F2F2F2', ec='#8A8A8A', fs=6.5):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.35,rounding_size=1.2', fc=c, ec=ec, lw=0.7))
    ax.text(x + w / 2, y + h / 2, t, ha='center', va='center', fontsize=fs, linespacing=1.3)
def arr(x1, x2, y=11): ax.add_patch(FancyArrowPatch((x1, y), (x2, y), arrowstyle='-|>', mutation_scale=8, lw=0.8, color='#555'))
box(0, 13, 15, 6, 'Query  $rdh54\\Delta$', c='#FBE3EE', ec=ACT); box(0, 3, 15, 6, 'Query  control', c='#EDEDED')
ax.text(17.6, 11, '×', fontsize=13, ha='center', va='center')
box(20.5, 6, 15, 10, 'Deletion library\n3,403 strains')
arr(36.5, 40.5); box(41.5, 6, 14, 10, 'Double mutants\n(SGA selection)'); arr(56.5, 60.5)
box(61.5, 13, 14, 6, 'YES'); box(61.5, 3, 14, 6, 'YES + 15 mM HU', c='#FFF3D6', ec=HU); arr(76.5, 80.5)
box(81.5, 5.5, 19, 11, 'Colony size, 72 h\n3 reps $rdh54\\Delta$, 2 reps ctrl\n\nscore = $rdh54\\Delta$ − ctrl', fs=5.9)
save(f, 'Fig4A_SGA_scheme')

# ---------------- B  YES vs HU interaction
E = G[G.evaluable_HU & (G['4601_YES_n'] >= 1) & (G['3206_YES_n'] >= 1)].copy()
q5 = G.loc[G.evaluable_HU, 'score_HU'].quantile(0.05)
f, ax = plt.subplots(figsize=(3.1, 3.0))
ax.axhspan(-5.2, q5, color='#FFF3D6', zorder=0, lw=0)
ax.scatter(E.score_YES, E.score_HU, s=4, color=GREY, lw=0, rasterized=True, zorder=1)
A = E[E['Actin cytoskeleton']]
ax.scatter(A.score_YES, A.score_HU, s=7, color=ACT, alpha=0.35, lw=0, zorder=2)
for g, (c, off) in HL.items():
    r = E[E.gene == g]
    if not len(r): continue
    r = r.iloc[0]; big = g == 'myo51'
    ax.scatter(r.score_YES, r.score_HU, s=30 if big else 18, color=c, edgecolor='k', lw=0.5, zorder=4)
    ax.annotate(it(g), (r.score_YES, r.score_HU), xytext=off, textcoords='offset points', fontsize=6.3, color=c,
                fontweight='bold' if big else 'normal', va='center',
                arrowprops=dict(arrowstyle='-', lw=0.4, color=c, shrinkA=0, shrinkB=2))
ax.axhline(0, lw=0.5, color='#999', zorder=1); ax.axvline(0, lw=0.5, color='#999', zorder=1)
ax.set_xlim(-8.8, 3.6); ax.set_ylim(-5.0, 2.8)
ax.text(3.5, q5 - 0.12, 'top 5% aggravating\nin HU', fontsize=5.5, color='#9A6B00', ha='right', va='top')
ax.set_xlabel('Interaction with $rdh54\\Delta$ in YES\n($log_2$ colony size, $rdh54\\Delta$ − control)')
ax.set_ylabel('Interaction with $rdh54\\Delta$ in YES + HU')
h = [plt.Line2D([], [], marker='o', ls='', color=c, mec='k', mew=0.4, ms=4, label=l) for c, l in
     [(ACT, 'actin cytoskeleton'), (MYO, 'type V myosin'), (REP, 'DNA repair'), (SPB, 'SPB / telomere')]]
ax.legend(handles=h, fontsize=5.5, frameon=False, loc='upper left', handletextpad=0.2, borderaxespad=0.2)
ax.set_title(f'n = {len(E):,} deletions', fontsize=6.5)
save(f, 'Fig4B_interaction_YES_vs_HU')
E.to_csv('output/SGA_screen/source_data_fig4B.csv')

# ---------------- C  ranked interaction in HU
H = G[G.evaluable_HU].sort_values('score_HU').reset_index(); H['rank'] = np.arange(1, len(H) + 1)
f, ax = plt.subplots(figsize=(3.5, 3.0))
ax.axvspan(0, len(H) * 0.05, color='#FFF3D6', lw=0, zorder=0)
ax.scatter(H['rank'], H.score_HU, s=2.5, color=GREY, lw=0, rasterized=True, zorder=1)
AH = H[H['Actin cytoskeleton']]
ax.scatter(AH['rank'], AH.score_HU, s=7, color=ACT, alpha=0.35, lw=0, zorder=2)
labs = [('rho4', ACT), ('msh2', REP), ('rdi1', ACT), ('ssm4', ACT), ('mcp6', SPB), ('bqt1', SPB), ('cap1', ACT),
        ('hus1', REP), ('myo51', MYO), ('mid1', ACT)]
ytxt = np.linspace(-4.7, -0.9, len(labs))
for (g, c), yt in zip(labs, ytxt):
    r = H[H.gene == g].iloc[0]; big = g == 'myo51'
    ax.scatter(r['rank'], r.score_HU, s=30 if big else 16, color=c, edgecolor='k', lw=0.5, zorder=4)
    ax.annotate(f"{it(g)}  #{int(r['rank'])}", (r['rank'], r.score_HU), xytext=(560, yt), textcoords='data',
                fontsize=6.2, color=c, va='center', fontweight='bold' if big else 'normal',
                arrowprops=dict(arrowstyle='-', lw=0.4, color=c, shrinkA=1, shrinkB=2))
ax.axhline(0, lw=0.5, color='#999')
ax.set_xlim(-20, len(H) + 20); ax.set_ylim(-5.0, 2.8)
ax.set_xlabel(f'Rank in YES + HU (of {len(H):,} deletions)'); ax.set_ylabel('Interaction with $rdh54\\Delta$ in YES + HU')
SL = G[G.absent_in_rdh54_HU & G['DNA repair'] & (G['4601_YES_n'] >= 1) & (G.gene != '0')].gene.tolist()
ax.text(0.98, 0.98, 'No colony in $rdh54\\Delta$ + HU\n(all replicates; present in control):\n' +
        '\n'.join(', '.join(it(g) for g in sorted(SL)[i:i+4]) for i in range(0, len(SL), 4)), transform=ax.transAxes, ha='right', va='top', fontsize=5.6, color=REP,
        bbox=dict(fc='white', ec=REP, lw=0.5, boxstyle='round,pad=0.3'))
save(f, 'Ranks_HU_not_in_figure')
H.to_csv('output/SGA_screen/source_data_ranks_HU.csv', index=False)

# ---------------- D  replicate values for actin-system genes
genes = ['rho4', 'rdi1', 'cap1', 'myo51']
f, axs = plt.subplots(1, len(genes), figsize=(5.2, 1.9), sharey=True); rows = []
groups = [('3206', 'YES'), ('4601', 'YES'), ('3206', 'YESHU'), ('4601', 'YESHU')]
xs = [0, 0.8, 2.0, 2.8]
for ax, g in zip(axs, genes):
    st = G[G.gene == g].index[0]
    x = d[d.strain == st].groupby(['query', 'medium', 'rep']).v.mean().reset_index()
    ax.axvspan(1.5, 3.3, color='#FFF3D6', lw=0, zorder=0)
    for xi, (q, m) in zip(xs, groups):
        v = x[(x['query'] == q) & (x.medium == m)].dropna(subset=['v'])
        c = (MYO if g == 'myo51' else ACT) if q == '4601' else '#8C8C8C'
        ax.scatter(np.full(len(v), xi) + np.linspace(-0.1, 0.1, len(v)) if len(v) > 1 else [xi] * len(v), v.v,
                   s=13, color=c, edgecolor='k', lw=0.3, zorder=3)
        if len(v): ax.plot([xi - 0.25, xi + 0.25], [v.v.mean()] * 2, color='k', lw=0.9, zorder=4)
        for _, r in v.iterrows(): rows.append(dict(gene=g, strain=st, query=q, medium=m, replicate=r.rep, log2_colony_size=r.v))
    ax.axhline(0, lw=0.5, color='#BBB', zorder=1)
    ax.set_xticks(xs); ax.set_xticklabels(['ctrl', '$rdh54\\Delta$'] * 2, rotation=55, ha='right', fontsize=5.8)
    ax.set_xlim(-0.5, 3.3); ax.set_title(it(g), color=MYO if g == 'myo51' else ACT, fontsize=7.5)
    ax.text(0.4, -6.6, 'YES', ha='center', fontsize=5.8, color='#666'); ax.text(2.4, -6.6, '+ HU', ha='center', fontsize=5.8, color='#9A6B00')
axs[0].set_ylim(-7.0, 1.2); axs[0].set_ylabel('$log_2$ colony size (normalized)')
for ax in axs[1:]: ax.tick_params(left=False)
f.subplots_adjust(wspace=0.12); save(f, 'Fig4C_replicates_actin')
pd.DataFrame(rows).to_csv('output/SGA_screen/source_data_fig4C.csv', index=False)
print('done; synthetic-lethal DNA-repair candidates:', SL)
