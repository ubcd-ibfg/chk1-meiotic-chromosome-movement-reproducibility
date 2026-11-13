"""Supp. Fig. 4 – reproducibility of the rdh54Δ Motif C phenotype across four independent series.
Each series: rdh54Δ vs the control imaged in the SAME batch (Motif C counts from J. Pinto-Cruz summaries).
Test: exact conditional test for two Poisson rates (binomial k1 | k1+k0, p = n1/(n1+n0)); exact 95% CI of the rate ratio.
Run:  python code/16_suppfig4_reproducibility.py  (from any directory)
"""
import os as _os; _os.chdir(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))  # run from the Reproducibility folder
import os, logging, numpy as np, pandas as pd, matplotlib
logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.stats import binomtest, beta
plt.rcParams.update({'font.family': ['Arial', 'Liberation Sans'], 'font.size': 7, 'axes.spines.top': False,
    'axes.spines.right': False, 'pdf.fonttype': 42, 'svg.fonttype': 'none', 'axes.linewidth': 0.7})
R = 'raw_data/S_pombe_imaging/submission_data_AFA_2026/'
SERIES = [  # label, figures where the series is shown, summary file, rdh54 row
    ('Series 1', 'Fig. 1; Fig. 2E', 'Figure_1/Figure_1_motifs_and_windows_summary.csv', 'rdh54'),
    ('Series 2', 'Fig. 2J–L; Fig. 3J, K', 'Figure_rdh54_double_mutants/Figure_rdh54_double_mutants_motifs_and_windows_summary.csv', 'rdh54'),
    ('Series 3', 'Fig. 4F, G', 'Figure_Rodrigo/Figure_Rodrigo_motifs_and_windows_summary.csv', 'rdh54_hyg'),
    ('Series 4', 'Fig. 4D, E', 'Figure_Microtubules_Actin/Figure_Microtubules_Actin_latA_motifs_and_windows_summary.csv', 'rdh54_nolatA')]

def rate_ratio(k1, n1, k0, n0):
    """Exact conditional test and CI for (k1/n1)/(k0/n0)."""
    n, p0 = k1 + k0, n1 / (n1 + n0)
    P = binomtest(int(k1), int(n), p0).pvalue
    lo_q = beta.ppf(0.025, k1, k0 + 1) if k1 > 0 else 0.0
    hi_q = beta.ppf(0.975, k1 + 1, k0) if k0 > 0 else 1.0
    conv = lambda q: np.inf if q >= 1 else q / (1 - q) * n0 / n1
    rr = np.inf if k0 == 0 else (k1 / n1) / (k0 / n0)
    return rr, conv(lo_q), conv(hi_q), P

rows = []
for lab, figs, f, s in SERIES:
    d = pd.read_csv(R + f)
    r = d[d.strain == s].iloc[0]; c = d[d.strain.str.lower().str.startswith('control')]
    n0 = int(c[c.motifC_mutant.isna()].N_cells.iloc[0])          # cells of the batch control
    w0 = float(c[c.motifC_mutant.isna()].windows.iloc[0])        # observation windows of the batch control
    k1, n1, k0, w1 = int(r.motifC_mutant), int(r.N_cells), int(r.motifC_control), float(r.windows)
    rr, lo, hi, P = rate_ratio(k1, n1, k0, n0)
    rows.append(dict(series=lab, shown_in=figs, source_file=f, rdh54_row=s, rdh54_events=k1, rdh54_cells=n1,
                     control_events=k0, control_cells=n0, rdh54_rate=k1 / n1, control_rate=k0 / n0,
                     rate_ratio=rr, CI_low=lo, CI_high=hi, P=P,
                     rdh54_windows=w1, control_windows=w0, rdh54_windows_per_cell=w1 / n1, control_windows_per_cell=w0 / n0,
                     rdh54_events_per_window=k1 / w1, control_events_per_window=k0 / w0))
T = pd.DataFrame(rows)
K1, N1, K0, N0 = T.rdh54_events.sum(), T.rdh54_cells.sum(), T.control_events.sum(), T.control_cells.sum()
W1, W0 = T.rdh54_windows.sum(), T.control_windows.sum()
rr, lo, hi, P = rate_ratio(K1, N1, K0, N0)
T = pd.concat([T, pd.DataFrame([dict(series='Pooled', shown_in='', source_file='sum of series 1–4', rdh54_row='',
        rdh54_events=K1, rdh54_cells=N1, control_events=K0, control_cells=N0, rdh54_rate=K1 / N1, control_rate=K0 / N0,
        rate_ratio=rr, CI_low=lo, CI_high=hi, P=P, rdh54_windows=W1, control_windows=W0, rdh54_windows_per_cell=W1 / N1,
        control_windows_per_cell=W0 / N0, rdh54_events_per_window=K1 / W1, control_events_per_window=K0 / W0)])], ignore_index=True)
os.makedirs('output/SuppFig4/panels', exist_ok=True); os.makedirs('Source_Data', exist_ok=True)
T.to_csv('output/SuppFig4/suppfig4_series.csv', index=False)
print(T[['series', 'rdh54_events', 'rdh54_cells', 'control_events', 'control_cells', 'rate_ratio', 'CI_low', 'CI_high', 'P']].round(4).to_string(index=False))
# Observation time (Discussion): rdh54Δ in the three series with ~1 event per cell (Series 1–3; Series 4 shares its motif
# dictionary with latrunculin A-treated cells, see Results). Events per cell vs events per observation window.
O = T.iloc[:3][['series', 'rdh54_events', 'rdh54_cells', 'rdh54_windows', 'rdh54_windows_per_cell', 'rdh54_rate', 'rdh54_events_per_window']].rename(
    columns={'rdh54_rate': 'rdh54_events_per_cell'})
O = pd.concat([O, pd.DataFrame([dict(series='max / min (Series 1–3)',
    rdh54_windows_per_cell=O.rdh54_windows_per_cell.max() / O.rdh54_windows_per_cell.min(),
    rdh54_events_per_cell=O.rdh54_events_per_cell.max() / O.rdh54_events_per_cell.min(),
    rdh54_events_per_window=O.rdh54_events_per_window.max() / O.rdh54_events_per_window.min())])], ignore_index=True)
print(O.round(3).to_string(index=False))

def fmtp(p): return 'P < 0.0001' if p < 1e-4 else f'P = {p:.3f}' if p >= 0.001 else f'P = {p:.4f}'
GREY, MAG = '#BDBDBD', '#E5407F'
f, (a, b) = plt.subplots(1, 2, figsize=(7.0, 2.4), gridspec_kw=dict(width_ratios=[1, 1.35], wspace=0.75))
# ---- A: Motif C events per cell, rdh54Δ vs same-batch control
S = T.iloc[:4]; x = np.arange(4); w = 0.36
a.bar(x - w / 2, S.control_rate, w, color=GREY, label='control (same batch)')
a.bar(x + w / 2, S.rdh54_rate, w, color=MAG, label='$rdh54\\Delta$')
for i, r in S.iterrows():
    a.text(i - w / 2, r.control_rate + 0.03, f'{r.control_events}/{r.control_cells}', ha='center', fontsize=5.5)
    a.text(i + w / 2, r.rdh54_rate + 0.03, f'{r.rdh54_events}/{r.rdh54_cells}', ha='center', fontsize=5.5)
    a.text(i, 1.62, fmtp(r.P).replace('P = ', ''), ha='center', fontsize=5.8, color='k' if r.P < 0.05 else '#888', rotation=0)
a.set_xticks(x); a.set_xticklabels(S.series); a.set_ylim(0, 1.72); a.set_yticks(np.arange(0, 1.21, 0.2)); a.spines['left'].set_bounds(0, 1.2); a.set_ylabel('Motif C events per cell')
a.legend(frameon=False, fontsize=6, loc='upper left', bbox_to_anchor=(0.0, 0.93), ncol=2, handlelength=1.2, columnspacing=1.0)
a.text(0.02, 0.76, 'numbers on bars: events / cells', transform=a.transAxes, ha='left', fontsize=5.3, color='#666')
a.text(-0.5, 1.62, 'P', fontsize=6, ha='right', style='italic')
# ---- B: forest plot of rate ratios (log scale); infinite bounds drawn as arrows
XMAX = 300; y = np.arange(5)[::-1]
b.axvline(1, ls='--', lw=0.6, color='#888')
for yi, (_, r) in zip(y, T.iterrows()):
    pooled = r.series == 'Pooled'; col = '#333' if pooled else MAG
    hi_ = min(r.CI_high, XMAX)
    b.plot([max(r.CI_low, 0.1), hi_], [yi, yi], color=col, lw=1.6 if pooled else 1.2, alpha=0.9, solid_capstyle='butt')
    if np.isinf(r.CI_high): b.annotate('', xy=(XMAX * 1.25, yi), xytext=(XMAX * 0.9, yi), arrowprops=dict(arrowstyle='-|>', color=col, lw=1.2, mutation_scale=7))
    if np.isinf(r.rate_ratio):
        b.text(max(r.CI_low, 0.1) * 1.08, yi + 0.2, '0 control events: ratio not finite', fontsize=5.3, color='#777', va='bottom')
    else:
        b.plot(r.rate_ratio, yi, marker='D' if pooled else 'o', ms=6 if pooled else 5, color=col, mec='white' if not pooled else col, mew=0.5, zorder=3)
    lab = '∞' if np.isinf(r.rate_ratio) else f'{r.rate_ratio:.1f}'
    ci = f'[{r.CI_low:.2f}–' + ('∞' if np.isinf(r.CI_high) else f'{r.CI_high:.1f}') + ']'
    b.text(XMAX * 2.2, yi, f'{lab} {ci}\n{fmtp(r.P)}', fontsize=5.8, va='center', fontweight='bold' if pooled else 'normal')
labels = [f'{r.series}\n{r.shown_in}' if r.series != 'Pooled' else f'Pooled\n{K1}/{N1} vs {K0}/{N0} cells' for _, r in T.iterrows()]
b.set_yticks(y); b.set_yticklabels(labels, fontsize=6.2); b.set_xscale('log'); b.set_xlim(0.1, XMAX * 1.4); b.set_ylim(-0.6, 4.6)
b.set_xticks([0.1, 1, 10, 100]); b.set_xticklabels(['0.1', '1', '10', '100'])
b.set_xlabel('Motif C rate ratio, $rdh54\\Delta$ / same-batch control (exact 95% CI)')
for ax, l in ((a, 'A'), (b, 'B')): ax.text(-0.32 if l == 'B' else -0.2, 1.06, l, transform=ax.transAxes, fontsize=9, fontweight='bold')
for e in ('pdf', 'svg'): f.savefig(f'output/SuppFig4/panels/SuppFig4_reproducibility.{e}', bbox_inches='tight')
f.savefig('output/SuppFig4/panels/SuppFig4_reproducibility.png', dpi=300, bbox_inches='tight')
with pd.ExcelWriter('Source_Data/Source_Data_SuppFig4.xlsx') as w_:
    pd.DataFrame({'Source data – Supplementary Fig. 4 (reproducibility of rdh54Δ Motif C)': [
        'Motif C counts from J. Pinto-Cruz *_motifs_and_windows_summary.csv files (raw_data/S_pombe_imaging/submission_data_AFA_2026/), unmodified.',
        'Each rdh54Δ series is compared with the control imaged in the same batch (control row N_cells; motifC_control of the rdh54 row).',
        'Exact conditional test for two Poisson rates (binomial k1 | k1+k0 with p = n1/(n1+n0)); exact 95% CI for the rate ratio (Clopper–Pearson on the conditional binomial).',
        'Series with 0 control events have a non-finite rate ratio; only the lower CI bound is informative.',
        'Pooled: sum of events, cells and windows across the four series.',
        'Observation window (Methods): number of 12-frame windows with t < −100 min (before meiosis I onset), summed over the cells of a strain; column "windows" of the summary files. Fractional values are as given in the summary files.',
        '*_windows_per_cell = windows / cells (observation time per cell); *_events_per_window = Motif C events / windows (events per unit of observation).',
        'Sheet Observation_time: rdh54Δ in Series 1–3, the three series compared in the Discussion (observation time per cell differs ×1.75, events per cell are ~1, events per window differ ~×2).',
        'Code: code/16_suppfig4_reproducibility.py']}).to_excel(w_, sheet_name='README', index=False)
    T.to_excel(w_, sheet_name='Series_and_pooled', index=False)
    O.to_excel(w_, sheet_name='Observation_time', index=False)
