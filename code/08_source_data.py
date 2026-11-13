"""Assemble Source_Data_S_pombe_imaging.xlsx (Fig. 1, 2E, J–L, 3J, K, 4D–G; Supp. Figs 1–4; Supp. Table 1) from output/ tables."""
from common import *
import glob
R = lambda f: pd.read_csv(os.path.join(OUT, f))
cells = R('cells_series1.csv'); stats = R('stats_series1_vs_control.csv')
mot = R('motif_rates_vs_same_batch_control.csv'); mot['q_BH_within_file_and_motif'] = mot.groupby(['summary_file', 'motif']).P.transform(lambda p: bh(p.fillna(1)))
b = R('suppfig2b_rf_vs_control.csv'); b['q_BH'] = bh(b.p_perm); b = b.merge(R('suppfig2b_ungrouped.csv').rename(columns={'genotype': 'g'}), on='g', how='left')
readme = pd.DataFrame({'Sheet': ['Fig1_cells', 'Fig1_SuppFig1_stats', 'Fig1f_zscores', 'Fig1g_PCA_loadings', 'SuppFig1a_cell_size', 'SuppFig1a_stats', 'SuppTable1',
    'SuppFig2a', 'SuppFig2b', 'SuppFig2c', 'SuppFig2d_pairwise', 'Motif_rates', 'Motif_key_comparisons', 'Fig2e', 'Fig2l'],
  'Content': ['Per-cell descriptors (Series 1), incl. PC1/PC2 (Fig. 1c, e, g; Supp. Fig. 1b–h)',
    'Mann–Whitney U vs control per variable and genotype; q = Benjamini–Hochberg across 13 genotypes (Fig. 1c, e; Supp. Fig. 1b–h)',
    'z-scores of per-cell means vs control (Fig. 1f)', 'PCA loadings and variance explained (Fig. 1g)',
    'Meiocyte size, one row per cell (Supp. Fig. 1a)', 'Supp. Fig. 1a statistics (BH across 13 genotypes)', 'Supplementary Table 1 (Series 1)',
    'Accuracy of nine classifiers, CV grouped by cell (Supp. Fig. 2a)', 'RF accuracy grouped by cell, permutation null mean, P (200 permutations), q (BH), ungrouped accuracy (grey ticks) (Supp. Fig. 2b)',
    'RF Gini importance (Supp. Fig. 2c)', 'Pairwise RF accuracy, 91 pairs (Supp. Fig. 2d)',
    'Motif A/C events and cells for every strain vs its same-batch control; rate ratio, exact 95% CI, exact conditional P, BH q; observation windows (12-frame windows with t < −100 min, summed over cells) for strain and control, windows per cell, events per window (Fig. 2e, j, k; Fig. 3j, k; Fig. 4d–g)',
    'Double mutants vs rdh54Δ and pooled rdh54Δ series (Fig. 2l; Supp. Fig. 4)', 'Volcano coordinates (Fig. 2e)', 'Motif C rate ratio vs rdh54Δ (Fig. 2l)']})
os.makedirs(os.path.join(ROOT, 'Source_Data'), exist_ok=True)
with pd.ExcelWriter(os.path.join(ROOT, 'Source_Data', 'Source_Data_S_pombe_imaging.xlsx')) as w:
    readme.to_excel(w, sheet_name='README', index=False)
    cells.to_excel(w, sheet_name='Fig1_cells', index=False); stats.to_excel(w, sheet_name='Fig1_SuppFig1_stats', index=False)
    R('fig1f_zscores.csv').to_excel(w, sheet_name='Fig1f_zscores', index=False); R('fig1g_pca_loadings.csv').to_excel(w, sheet_name='Fig1g_PCA_loadings', index=False)
    R('suppfig1a_cell_size.csv').to_excel(w, sheet_name='SuppFig1a_cell_size', index=False); R('stats_suppfig1a_cell_size.csv').to_excel(w, sheet_name='SuppFig1a_stats', index=False)
    R('supp_table1_series1.csv').to_excel(w, sheet_name='SuppTable1', index=False)
    R('suppfig2a_benchmark.csv').to_excel(w, sheet_name='SuppFig2a', index=False); b.to_excel(w, sheet_name='SuppFig2b', index=False)
    R('suppfig2c_gini.csv').to_excel(w, sheet_name='SuppFig2c', index=False); R('suppfig3_pairwise.csv').to_excel(w, sheet_name='SuppFig2d_pairwise', index=False)
    mot.to_excel(w, sheet_name='Motif_rates', index=False); R('motif_key_comparisons.csv').to_excel(w, sheet_name='Motif_key_comparisons', index=False)
    pd.read_csv(os.path.join(OUT, 'motif_panels', 'Fig2e_data.csv')).to_excel(w, sheet_name='Fig2e', index=False)
    pd.read_csv(os.path.join(OUT, 'motif_panels', 'Fig2l_data.csv')).to_excel(w, sheet_name='Fig2l', index=False)
print('ok')
