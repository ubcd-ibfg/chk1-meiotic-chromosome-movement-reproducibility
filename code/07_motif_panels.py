"""Motif-frequency panels: Fig. 2e (volcano), Fig. 2j–l (rdh54Δ double mutants), Fig. 3j, k (microtubule series),
Fig. 4d, e (latrunculin A series), Fig. 4f, g (myo51Δ series). Run after 03_motif_statistics.py.
grey = control decomposed in the same batch; magenta = strain; error bars, exact Poisson 95 % CI;
P, exact conditional test for two Poisson rates (cells as unit of exposure).
Fig. 3k also shows the bracket rdh54Δ hrs1-T16A vs rdh54Δ (Motif C), as in the paper."""
from common import *
import logging; logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)
from scipy.stats import binomtest, chi2, beta
from matplotlib.lines import Line2D
plt.rcParams.update({'font.size': 7, 'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5})
T = pd.read_csv(os.path.join(OUT, 'motif_rates_vs_same_batch_control.csv'))
V = os.path.join(OUT, 'motif_panels'); os.makedirs(V, exist_ok=True)
CTRL, MUT = '#B0B0B0', '#D6338A'; CA, CC_ = '#E07B73', '#4FAE3E'
CATC = {'DSB formation': '#E8743B', 'DNA repair': '#D6338A', 'Checkpoint': '#6A5ACD', 'Chromatin': '#2A9D6F', 'Cohesion': '#D9A21B', 'Treatment': '#7F7F7F'}
CATG = {'rec12': 'DSB formation', 'rdh54': 'DNA repair', 'exo1': 'DNA repair', 'rev3': 'DNA repair', 'fml1': 'DNA repair', 'tel1': 'Checkpoint', 'rad3': 'Checkpoint',
        'cds1': 'Checkpoint', 'chk1': 'Checkpoint', 'mek1': 'Checkpoint', 'dsk1': 'Checkpoint', 'snf22': 'Chromatin', 'rec8': 'Cohesion'}
def pc(k1, n1, k0, n0):
    k1, k0 = int(k1), int(k0); return np.nan if k1 + k0 == 0 else binomtest(k1, k1 + k0, n1 / (n1 + n0)).pvalue
def poisci(k, n): return (0 if k == 0 else chi2.ppf(0.025, 2 * k) / 2) / n, chi2.ppf(0.975, 2 * k + 2) / 2 / n
def rrci(k1, n1, k0, n0):
    lo = beta.ppf(0.025, k1, k0 + 1) if k1 > 0 else 0; hi = beta.ppf(0.975, k1 + 1, k0) if k0 > 0 else 1
    f = lambda q: q / (1 - q) * n0 / n1 if q < 1 else np.inf; return f(lo), f(hi)
def fp(p): return 'n.a.' if np.isnan(p) else 'P<0.0001' if p < 1e-4 else f'P={p:.3f}' if p >= 0.001 else f'P={p:.4f}'
def row(f, s, M): return T[(T.summary_file.str.startswith(f)) & (T.strain == s) & (T.motif == M)].iloc[0]
def vsave(fig, n):
    for e in ('pdf', 'svg'): fig.savefig(os.path.join(V, f'{n}.{e}'), bbox_inches='tight')
    plt.close(fig)
def paired(ax, f, keys, labels, M):
    x = 0; ticks = []; ymax = 0
    for k, lab in zip(keys, labels):
        r = row(f, k, M)
        for j, (kk, nn, col) in enumerate(((r.control_events_same_batch, r.control_cells, CTRL), (r.events, r.cells, MUT))):
            v = kk / nn; lo, hi = poisci(kk, nn); xx = x + j * 0.42
            ax.bar(xx, v, 0.38, color=col, zorder=2)
            ax.errorbar(xx, v, yerr=[[v - lo], [hi - v]], fmt='none', ecolor='#333', elinewidth=0.6, capsize=1.5, zorder=3); ymax = max(ymax, hi)
        p = r.P
        ax.text(x + 0.21, max(poisci(r.control_events_same_batch, r.control_cells)[1], poisci(r.events, r.cells)[1]) + 0.04, fp(p),
                ha='center', va='bottom', fontsize=6, color='k' if p < 0.05 else '#999', rotation=45)
        ticks.append((x + 0.21, lab)); x += 1.15
    ax.set_xticks([t for t, _ in ticks]); ax.set_xticklabels([l for _, l in ticks], rotation=40, ha='right', fontstyle='italic')
    ax.set_ylabel(f'Motif {M} events per cell'); ax.set_title(f'Motif {M}', color=CA if M == 'A' else CC_); ax.set_ylim(0, ymax * 1.35)
    return ax
# ---- Fig. 2e volcano
fig, axs = plt.subplots(1, 2, figsize=(5.4, 2.4)); rows = []
sub = T[T.summary_file.str.startswith(('Figure_1_', 'Supplementary_Figure_1'))]
for ax, M in zip(axs, 'AC'):
    pts = {}
    for _, r in sub[sub.motif == M].iterrows():
        tr = r.summary_file.startswith('Supplementary'); p = pc(r.events, r.cells, r.control_events_same_batch, r.control_cells)
        lr = np.log2(((r.events + 0.5) / r.cells) / ((r.control_events_same_batch + 0.5) / r.control_cells)); yl = -np.log10(p) if not np.isnan(p) else 0
        cat = 'Treatment' if tr else CATG[r.strain]; sig = p < 0.05
        lab = {'ZEO': 'Zeocin'}.get(r.strain, r.strain if tr else r.strain + 'Δ')
        ax.scatter(lr, yl, s=24 if sig else 14, marker='D' if tr else 'o', color=CATC[cat] if sig else '#CFCFCF', edgecolor='white', lw=0.5, zorder=3)
        rows.append(dict(motif=M, label=lab, log2_rate_ratio=lr, minus_log10_P=yl, P=p))
        pts.setdefault((round(lr, 2), round(yl, 2)), []).append((lab, cat, sig, tr))
    for (lr, yl), L in pts.items():
        if yl < 0.75 and not L[0][2]: continue
        txt = ', '.join(sorted(l for l, _, _, _ in L)); cat, sig = L[0][1], L[0][2]
        off = {'snf22Δ': (-30, -9), 'Zeocin': (-16, 5), 'HU, rev3Δ': (-14, 5), 'MMS': (4, -8), 'rad3Δ, snf22Δ, tel1Δ': (-24, 6), 'mek1Δ': (4, 2)}.get(txt, (4, 3))
        ax.annotate(txt, (lr, yl), xytext=off, textcoords='offset points', fontsize=6, fontstyle='italic', color=CATC[cat] if sig else '#777')
    ax.axhline(-np.log10(0.05), ls='--', lw=0.6, color='#999'); ax.axvline(0, lw=0.6, color='#666')
    ax.set_xlabel(f'log$_2$ Motif {M} rate vs same-batch control'); ax.set_ylabel('−log$_{10}$ P (exact conditional)')
    ax.set_title(f'Motif {M}', color=CA if M == 'A' else CC_)
fig.tight_layout(); vsave(fig, 'Fig2e_volcano'); pd.DataFrame(rows).to_csv(os.path.join(V, 'Fig2e_data.csv'), index=False)
# ---- Fig. 2j–l double mutants
f = 'Figure_rdh54_double'; keys = ['rdh54', 'rdh54_rec12', 'rdh54_cds1', 'rdh54_rad3', 'rdh54_chk1']
labs = ['rdh54Δ', 'rdh54Δ rec12Δ', 'rdh54Δ cds1Δ', 'rdh54Δ rad3Δ', 'rdh54Δ chk1Δ']
for M, n in (('A', 'Fig2j_doubles_motifA'), ('C', 'Fig2k_doubles_motifC')):
    fig, ax = plt.subplots(figsize=(2.3, 2.3)); paired(ax, f, keys, labs, M); vsave(fig, n)
fig, ax = plt.subplots(figsize=(2.6, 2.0)); r0 = row(f, 'rdh54', 'C'); rr_rows = []
for i, (k, l) in enumerate(zip(keys[1:], labs[1:])):
    r = row(f, k, 'C'); rr = (r.events / r.cells) / (r0.events / r0.cells); lo, hi = rrci(r.events, r.cells, r0.events, r0.cells); p = pc(r.events, r.cells, r0.events, r0.cells)
    ax.plot([max(lo, 0.003), hi], [i, i], color='#555', lw=1); ax.scatter(rr, i, s=22, color=MUT if p < 0.05 else '#AAA', zorder=3)
    ax.text(3.3, i, fp(p), va='center', fontsize=6, color='k' if p < 0.05 else '#999'); rr_rows.append(dict(comparison=f'{l} vs rdh54Δ', rate_ratio=rr, CI_low=lo, CI_high=hi, P=p))
ax.axvline(1, ls='--', lw=0.6, color='#999'); ax.set_xscale('log'); ax.set_xlim(0.003, 8)
ax.set_yticks(range(4)); ax.set_yticklabels(labs[1:], fontstyle='italic'); ax.invert_yaxis()
ax.set_xlabel('Motif C rate ratio vs rdh54Δ (95% CI)'); ax.set_title('Double mutant vs rdh54Δ'); vsave(fig, 'Fig2l_rate_ratio')
pd.DataFrame(rr_rows).to_csv(os.path.join(V, 'Fig2l_data.csv'), index=False)
# ---- Fig. 3j, k / 4d, e / 4f, g
for f, keys, labs, names in [
    ('Figure_Microtubules_Actin_microtubules', ['rdh54', 'klp6', 'klp6_rdh54', 'hrs1_t16A', 'hrs1_t16A_rdh54'], ['rdh54Δ', 'klp6Δ', 'rdh54Δ klp6Δ', 'hrs1-T16A', 'rdh54Δ hrs1-T16A'], ('Fig3j_MT_motifA', 'Fig3k_MT_motifC')),
    ('Figure_Microtubules_Actin_latA', ['Control_latA', 'rdh54_nolatA', 'rdh54_latA'], ['control + latA', 'rdh54Δ + DMSO', 'rdh54Δ + latA'], ('Fig4d_latA_motifA', 'Fig4e_latA_motifC')),
    ('Figure_Rodrigo', ['rdh54_hyg', 'myo51_kan', 'myo51_kan_rdh54_hyg'], ['rdh54Δ', 'myo51Δ', 'rdh54Δ myo51Δ'], ('Fig4f_myo51_motifA', 'Fig4g_myo51_motifC'))]:
    for M, n in zip('AC', names):
        fig, ax = plt.subplots(figsize=(0.55 * len(keys) + 0.8, 2.3)); paired(ax, f, keys, labs, M)
        if n == 'Fig3k_MT_motifC':   # bracket: rdh54Δ hrs1-T16A vs rdh54Δ (Motif C, exact conditional test; same as motif_key_comparisons.csv)
            i1, i2 = keys.index('rdh54'), keys.index('hrs1_t16A_rdh54'); a, b = row(f, 'rdh54', M), row(f, 'hrs1_t16A_rdh54', M)
            x1, x2 = i1 * 1.15 + 0.42, i2 * 1.15 + 0.42                       # magenta (strain) bars
            y1 = poisci(a.events, a.cells)[1] + 0.03                            # top of rdh54Δ CI
            y2 = max(poisci(b.events, b.cells)[1], poisci(b.control_events_same_batch, b.control_cells)[1]) + 0.04 + 0.45   # above the rotated P label
            yb = max(y1, y2) + 0.25; p = pc(b.events, b.cells, a.events, a.cells)
            ax.plot([x1, x1, x2, x2], [y1, yb, yb, y2], color='k', lw=0.6, zorder=4)
            ax.text((x1 + x2) / 2, yb + 0.03, fp(p), ha='center', va='bottom', fontsize=6, color='k' if p < 0.05 else '#999')
            ax.set_ylim(0, max(ax.get_ylim()[1], yb * 1.15))
        vsave(fig, n)
print('saved to', V)
