"""Vector re-draw of Fig. 1G (PCA) and Supp. Fig. 2A–D at final print size (text ≥ 6 pt, no raster layers).
Inputs (written by 01/02): output/cells_series1.csv, fig1g_pca_loadings.csv, suppfig2a_benchmark.csv,
suppfig2b_rf_vs_control.csv, suppfig2c_gini.csv, suppfig3_pairwise.csv.   Outputs: output/vector_panels/*.pdf|svg. Also needs output/suppfig2b_ungrouped.csv (06).
"""
import os, sys, logging, numpy as np, pandas as pd, matplotlib
logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import ORDER, CAT
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
IN = os.environ.get('VEC_IN', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'output', 'S_pombe_imaging'))
OUT = os.path.join(IN, 'vector_panels'); os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({'font.family': ['Arial', 'Liberation Sans'], 'font.size': 7, 'axes.spines.top': False, 'axes.spines.right': False,
    'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'pdf.fonttype': 42, 'svg.fonttype': 'none'})
MIN = 6.0
def save(f, n):
    for e in ('pdf', 'svg'): f.savefig(os.path.join(OUT, f'{n}.{e}'), bbox_inches='tight')
    plt.close(f)
def ital(ax, axis='y'):
    for t in (ax.get_yticklabels() if axis == 'y' else ax.get_xticklabels()): t.set_fontstyle('italic')

# ---------------- Fig. 1G  PCA of individual cells
cell = pd.read_csv(os.path.join(IN, 'cells_series1.csv'))
D6 = ['area', 'major', 'minor', 'circularity', 'convexity', 'cal.speed']
X = StandardScaler().fit_transform(cell[D6]); pca = PCA().fit(X); S = pca.transform(X)
L = pd.read_csv(os.path.join(IN, 'fig1g_pca_loadings.csv'), index_col=0)
assert np.allclose(np.abs(pca.components_[:2].T), np.abs(L[['PC1', 'PC2']].values), atol=1e-6)
f, ax = plt.subplots(figsize=(2.75, 2.55))
for g in ORDER:
    m = (cell.genotype == g).values; ax.scatter(S[m, 0], S[m, 1], s=5, color=CAT[g], alpha=0.35, lw=0)
LB = ['Area', 'Major axis', 'Minor axis', 'Circularity', 'Convexity', 'Velocity']
POS = {'Area': (-2.3, 2.9), 'Minor axis': (1.6, 2.9), 'Major axis': (-4.2, 1.4), 'Circularity': (4.0, 1.6), 'Convexity': (4.5, 0.15), 'Velocity': (-4.4, -0.35)}
for j, l in enumerate(LB):
    x, y = pca.components_[0, j] * 2.8, pca.components_[1, j] * 2.8
    ax.annotate('', xy=(x, y), xytext=(0, 0), arrowprops=dict(arrowstyle='-|>', color='#7A7A7A', lw=0.7, mutation_scale=5), zorder=6)
    tx, ty = POS[l]
    ax.plot([x, tx], [y, ty], color='#AFAFAF', lw=0.4, ls=':', zorder=2)
    ax.text(tx, ty, l, fontsize=MIN, color='#555', ha='center', va='center', zorder=8, bbox=dict(fc='white', ec='none', alpha=0.85, pad=0.4))
for g in ORDER:
    m = (cell.genotype == g).values; c = S[m, :2].mean(0)
    ax.scatter(*c, s=30, color=CAT[g], edgecolor='k', lw=0.6, zorder=5)
lab = {'rec8Δ': (6, 0), 'mek1Δ': (2, 7), 'rdh54Δ': (6, -7), 'control': (-12, -17)}
for g, o in lab.items():
    c = S[(cell.genotype == g).values, :2].mean(0)
    ax.annotate(g, c, xytext=o, textcoords='offset points', fontsize=6.5, fontstyle='normal' if g == 'control' else 'italic', zorder=7)
ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0] * 100:.1f} %)'); ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1] * 100:.1f} %)')
ax.set_title('PCA of meiotic chromosome dynamics\n(individual cells)', fontsize=7)
save(f, 'Fig1G_PCA')

# ---------------- Supp. Fig. 2
GEN = ['rec8Δ', 'chk1Δ', 'mek1Δ', 'snf22Δ', 'rdh54Δ', 'rec12Δ', 'tel1Δ', 'exo1Δ', 'dsk1Δ', 'fml1Δ', 'rev3Δ', 'rad3Δ', 'cds1Δ']
CL = ['Logistic regression', 'Decision tree', 'Naive Bayes', 'k-NN (k=5)', 'SVM (RBF)', 'Random forest', 'LDA', 'QDA', 'MLP']
mag = plt.get_cmap('magma_r')
def heat(ax, M, cmap, vmin, vmax, fmt='{:.2f}', fs=MIN):
    h, w = M.shape
    pc = ax.pcolormesh(np.arange(w + 1) - 0.5, np.arange(h + 1) - 0.5, M, cmap=cmap, vmin=vmin, vmax=vmax, edgecolors='white', linewidth=0.3)
    for i in range(h):
        for j in range(w):
            v = M[i, j]
            if np.isnan(v): continue
            r = (v - vmin) / (vmax - vmin); ax.text(j, i, fmt.format(v), ha='center', va='center', fontsize=fs, color='white' if r > 0.6 else 'black')
    ax.set_xlim(-0.5, w - 0.5); ax.set_ylim(h - 0.5, -0.5)
    for s in ax.spines.values(): s.set_visible(False)
    ax.tick_params(length=0)
    return pc
# A – benchmark
A = pd.read_csv(os.path.join(IN, 'suppfig2a_benchmark.csv'), index_col=0).loc[GEN, CL]
f, ax = plt.subplots(figsize=(3.3, 2.75))
pc = heat(ax, A.values, mag, 0.4, 0.9)
ax.set_xticks(range(len(CL))); ax.set_xticklabels(CL, rotation=50, ha='right'); ax.set_yticks(range(len(GEN))); ax.set_yticklabels(GEN); ital(ax)
cb = f.colorbar(pc, ax=ax, shrink=0.6, pad=0.02); cb.set_label('Accuracy'); cb.outline.set_linewidth(0.5); cb.solids.set_rasterized(False)
ax.set_title('Nine classifiers, CV grouped by cell', fontsize=7)
save(f, 'SuppFig2A_classifiers'); A.to_csv(os.path.join(OUT, 'SuppFig2A_data.csv'))
# B – RF vs control with permutation P
B = pd.read_csv(os.path.join(IN, 'suppfig2b_rf_vs_control.csv')).set_index('g').loc[GEN]
f, ax = plt.subplots(figsize=(2.2, 2.75)); y = np.arange(len(GEN))
sig = B.p_perm < 0.05
U = pd.read_csv(os.path.join(IN, 'suppfig2b_ungrouped.csv')).set_index('genotype').loc[GEN]; B['acc_ungrouped'] = U.acc_ungrouped
ax.barh(y, B.acc_grouped - 0.5, left=0.5, color=np.where(sig, '#E5407F', '#D3D3D3'), height=0.62)
ax.scatter(B.acc_ungrouped, y, marker='|', s=55, color='#7F7F7F', lw=1.2, zorder=3)
for yi, (g, r) in zip(y, B.iterrows()):
    ax.text(max(r.acc_grouped, r.acc_ungrouped, 0.5) + 0.015, yi, f'P={r.p_perm:.3f}', va='center', fontsize=MIN, color='k' if r.p_perm < 0.05 else '#888', fontweight='bold' if r.p_perm < 0.05 else 'normal')
ax.axvline(0.5, ls='--', lw=0.6, color='k'); ax.set_yticks(y); ax.set_yticklabels(GEN); ital(ax); ax.set_ylim(len(GEN) - 0.5, -0.5)
ax.set_xlim(0.4, 1.0); ax.set_xlabel('Random forest accuracy'); ax.set_title('RF vs control (permutation P)', fontsize=7)
save(f, 'SuppFig2B_rf_vs_control'); B.to_csv(os.path.join(OUT, 'SuppFig2B_data.csv'))
# C – Gini importance
C = pd.read_csv(os.path.join(IN, 'suppfig2c_gini.csv')).set_index('g').loc[GEN, ['area', 'minor', 'circularity', 'cal.speed', 'ang.speed']]
f, ax = plt.subplots(figsize=(1.95, 2.75))
pc = heat(ax, C.values, plt.get_cmap('Blues'), 0.0, 0.4)
ax.set_xticks(range(5)); ax.set_xticklabels(['Area', 'Minor axis', 'Circularity', 'Velocity', 'Angular speed'], rotation=50, ha='right')
ax.set_yticks(range(len(GEN))); ax.set_yticklabels(GEN); ital(ax); ax.set_title('RF Gini importance', fontsize=7)
save(f, 'SuppFig2C_gini'); C.to_csv(os.path.join(OUT, 'SuppFig2C_data.csv'))
# D – pairwise
P = pd.read_csv(os.path.join(IN, 'suppfig3_pairwise.csv'))
M = pd.DataFrame(np.nan, index=ORDER, columns=ORDER)
for _, r in P.iterrows(): M.loc[r.a, r.b] = r.acc; M.loc[r.b, r.a] = r.acc
f, ax = plt.subplots(figsize=(4.3, 3.9))
pc = heat(ax, M.values, mag, 0.4, 0.9)
ax.set_xticks(range(len(ORDER))); ax.set_xticklabels(ORDER, rotation=50, ha='right'); ax.set_yticks(range(len(ORDER))); ax.set_yticklabels(ORDER)
ital(ax); ital(ax, 'x')
for t in ax.get_xticklabels() + ax.get_yticklabels():
    if t.get_text() == 'control': t.set_fontstyle('normal')
cb = f.colorbar(pc, ax=ax, shrink=0.55, pad=0.02); cb.set_label('RF accuracy (CV grouped by cell)'); cb.outline.set_linewidth(0.5); cb.solids.set_rasterized(False)
GR = [('DSB\nformation', 1, 1, '#E8743B'), ('DNA\nrepair', 2, 5, '#D6338A'), ('Checkpoint', 6, 11, '#6A5ACD'), ('Chromatin', 12, 12, '#2A9D6F'), ('Cohesion', 13, 13, '#D9A21B')]
for lab, a, b, c in GR:
    ax.plot([-2.75, -2.75], [a - 0.42, b + 0.42], color=c, lw=3, clip_on=False, solid_capstyle='butt')
    ax.text(-3.0, (a + b) / 2, lab, ha='right', va='center', fontsize=6.5, clip_on=False)
ax.set_title('Pairwise classification', fontsize=7)
save(f, 'SuppFig2D_pairwise'); M.to_csv(os.path.join(OUT, 'SuppFig2D_data.csv'))
print('saved to', OUT)
