"""Supp. Fig. 2 as one vector figure at final print size (174 mm wide); numbers 5 pt, labels ≥ 6 pt."""
from common import *
import logging; logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)
plt.rcParams.update({'font.family': ['Arial', 'Liberation Sans'], 'font.size': 6.5, 'axes.linewidth': 0.5,
                     'xtick.major.width': 0.5, 'ytick.major.width': 0.5, 'xtick.major.size': 2, 'ytick.major.size': 2})
MM = 1 / 25.4
GEN = ['rec8Δ', 'chk1Δ', 'mek1Δ', 'snf22Δ', 'rdh54Δ', 'rec12Δ', 'tel1Δ', 'exo1Δ', 'dsk1Δ', 'fml1Δ', 'rev3Δ', 'rad3Δ', 'cds1Δ']
CL = ['Logistic regression', 'Decision tree', 'Naive Bayes', 'k-NN (k=5)', 'SVM (RBF)', 'Random forest', 'LDA', 'QDA', 'MLP']
R = lambda f: pd.read_csv(os.path.join(OUT, f))
A = R('suppfig2a_benchmark.csv'); A = A.set_index(A.columns[0]).loc[GEN, CL]
B = R('suppfig2b_rf_vs_control.csv'); B = B.set_index(B.columns[0]).loc[GEN]
U = R('suppfig2b_ungrouped.csv').set_index('genotype').loc[GEN]
C = R('suppfig2c_gini.csv'); C = C.set_index(C.columns[0]).loc[GEN, ['area', 'minor', 'circularity', 'cal.speed', 'ang.speed']]
P = R('suppfig3_pairwise.csv'); P.columns = ['a', 'b', 'acc']
M = pd.DataFrame(np.nan, index=ORDER, columns=ORDER)
for _, r in P.iterrows(): M.loc[r.a, r.b] = r.acc; M.loc[r.b, r.a] = r.acc
mag = plt.get_cmap('magma_r')
def heat(ax, Mv, cmap, vmin, vmax, fs=5):
    h, w = Mv.shape
    pc = ax.pcolormesh(np.arange(w + 1) - 0.5, np.arange(h + 1) - 0.5, Mv, cmap=cmap, vmin=vmin, vmax=vmax, edgecolors='white', linewidth=0.25)
    for i in range(h):
        for j in range(w):
            v = Mv[i, j]
            if np.isnan(v): continue
            ax.text(j, i, f'{v:.2f}', ha='center', va='center', fontsize=fs, color='white' if (v - vmin) / (vmax - vmin) > 0.6 else 'black')
    ax.set_xlim(-0.5, w - 0.5); ax.set_ylim(h - 0.5, -0.5)
    for s in ax.spines.values(): s.set_visible(False)
    ax.tick_params(length=0)
    return pc
def ital(labels): return labels
fig = plt.figure(figsize=(174 * MM, 165 * MM))
# panel A
axA = fig.add_axes([0.085, 0.60, 0.26, 0.33])
pcA = heat(axA, A.values, mag, 0.4, 0.9)
axA.set_xticks(range(9)); axA.set_xticklabels(CL, rotation=55, ha='right', fontsize=6); axA.set_yticks(range(13)); axA.set_yticklabels(GEN, fontstyle='italic', fontsize=6.5)
axA.set_title('Nine classifiers, CV grouped by cell', fontsize=7)
cA = fig.add_axes([0.352, 0.68, 0.008, 0.18]); cb = fig.colorbar(pcA, cax=cA); cb.set_label('Accuracy', fontsize=6); cb.ax.tick_params(labelsize=5.5, length=1.5); cb.outline.set_linewidth(0.4); cb.solids.set_rasterized(False); cb.solids.set_edgecolor('face')
# panel B
axB = fig.add_axes([0.47, 0.60, 0.20, 0.33]); y = np.arange(13); sig = B.iloc[:, 2] < 0.05 if 'p_perm' not in B else B.p_perm < 0.05
acc = B.iloc[:, 0]; p = B.p_perm
axB.barh(y, acc - 0.5, left=0.5, color=np.where(p < 0.05, '#E5407F', '#D3D3D3'), height=0.62)
axB.scatter(U.acc_ungrouped, y, marker='|', s=22, color='#7F7F7F', lw=1.0, zorder=3)
for yi, a_, u_, pp in zip(y, acc, U.acc_ungrouped, p):
    axB.text(max(a_, u_, 0.5) + 0.012, yi, f'P={pp:.3f}', va='center', fontsize=5, color='k' if pp < 0.05 else '#888', fontweight='bold' if pp < 0.05 else 'normal')
axB.axvline(0.5, ls='--', lw=0.5, color='k'); axB.set_yticks(y); axB.set_yticklabels(GEN, fontstyle='italic', fontsize=6.5); axB.set_ylim(12.5, -0.5)
axB.set_xlim(0.4, 1.0); axB.set_xticks([0.4, 0.6, 0.8, 1.0]); axB.set_xlabel('Random forest accuracy', fontsize=6.5)
axB.spines['top'].set_visible(False); axB.spines['right'].set_visible(False); axB.set_title('RF vs control (permutation P)', fontsize=7)
# panel C
axC = fig.add_axes([0.80, 0.60, 0.17, 0.33])
heat(axC, C.values, plt.get_cmap('Blues'), 0.0, 0.4)
axC.set_xticks(range(5)); axC.set_xticklabels(['Area', 'Minor axis', 'Circularity', 'Velocity', 'Angular speed'], rotation=55, ha='right', fontsize=6)
axC.set_yticks(range(13)); axC.set_yticklabels(GEN, fontstyle='italic', fontsize=6.5); axC.set_title('RF Gini importance', fontsize=7)
# panel D
axD = fig.add_axes([0.30, 0.07, 0.42, 0.40])
pcD = heat(axD, M.values, mag, 0.4, 0.9)
axD.set_xticks(range(14)); axD.set_xticklabels(ORDER, rotation=55, ha='right', fontsize=6); axD.set_yticks(range(14)); axD.set_yticklabels(ORDER, fontsize=6.5)
for t in axD.get_xticklabels() + axD.get_yticklabels(): t.set_fontstyle('normal' if t.get_text() == 'control' else 'italic')
for lab, a, b, c in [('DSB\nformation', 1, 1, '#E8743B'), ('DNA\nrepair', 2, 5, '#D6338A'), ('Checkpoint', 6, 11, '#6A5ACD'), ('Chromatin', 12, 12, '#2A9D6F'), ('Cohesion', 13, 13, '#D9A21B')]:
    axD.plot([-2.6, -2.6], [a - 0.42, b + 0.42], color=c, lw=2.5, clip_on=False, solid_capstyle='butt')
    axD.text(-2.9, (a + b) / 2, lab, ha='right', va='center', fontsize=6.5, clip_on=False)
axD.set_title('Pairwise classification', fontsize=7)
cD = fig.add_axes([0.735, 0.17, 0.010, 0.20]); cb = fig.colorbar(pcD, cax=cD); cb.set_label('RF accuracy (CV grouped by cell)', fontsize=6); cb.ax.tick_params(labelsize=5.5, length=1.5); cb.outline.set_linewidth(0.4); cb.solids.set_rasterized(False); cb.solids.set_edgecolor('face')
for x, yy, l in [(0.005, 0.965, 'A'), (0.385, 0.965, 'B'), (0.70, 0.965, 'C'), (0.12, 0.50, 'D')]: fig.text(x, yy, l, fontsize=9, fontweight='bold')
os.makedirs(os.path.join(OUT, 'vector_panels'), exist_ok=True)
for e in ('pdf', 'svg'): fig.savefig(os.path.join(OUT, 'vector_panels', f'SuppFig2_composite_174mm.{e}'))
print('ok')
