"""Fig. 4 SGA – step 3: assemble Source_Data/Source_Data_Fig4_SGA.xlsx from output/SGA_screen/ (run after 11 and 12)."""
import os; os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))  # run from the Reproducibility folder
import pandas as pd
O = 'output/SGA_screen/'
README = [('Raw data', 'raw_data/data_heatmap_completo_replica123_sergio_log2fold.csv (label;strain;gene2;log2fold_value_normalized2), unmodified'), ('Query strains', '4601 = rdh54Δ query; 3206 = control query'), ('Media', 'YES and YES + 15 mM hydroxyurea (label "YESHU")'), ('Time point', '72 h (24 h and 48 h are in the raw file but not used)'), ('Replicates', 'rdh54Δ query: A, B, C; control query: A, C. Each replicate value = mean of duplicate spots'), ('Colony value', 'log2fold_value_normalized2 (colony size, plate-normalized log2)'), ('NA', 'NA = no colony detected (synthetic lethality or replicate/pinning failure); kept as missing, NOT set to 0'), ('Interaction score', 'score = mean(rdh54Δ replicates) − mean(control replicates), per medium; negative = rdh54Δ aggravates growth of the deletion'), ('Evaluable (HU)', '≥2 rdh54Δ replicates and ≥1 control replicate with a colony in YES+HU (1,747 deletions)'), ('Rank', 'rank_HU: 1 = most negative score_HU among evaluable deletions'), ('GO categories', 'PomBase GAF (2026-05-21) propagated through GO is_a/part_of (go-basic 2026-07-26). Actin: GO:0030036, 0015629, 0003779, 0007015; Myosin: GO:0016459; DNA repair: GO:0006281; Microtubule: GO:0000226, 0015630, 0008017; DNA replication: GO:0006260'), ('Absent in rdh54 HU', 'no colony in any rdh54Δ YES+HU replicate but colony in control YES+HU'), ('Original ranking', 'Top_200_SGA.R (NA→0, rdh54Δ HU 72 h only, top-500 in ≥2 replicates, sum/3). Reproduced: 199/200 identical to genes_4601_200_new.txt'), ('Code', 'code/11_sga_scores.py, 12_sga_panels.py, 17_sga_source_data.py')]
G = pd.read_csv(O + 'sga_gene_scores.csv')
sheets = {
    'Fig4B_YES_vs_HU': pd.read_csv(O + 'source_data_fig4B.csv'),
    'Ranks_HU': pd.read_csv(O + 'source_data_ranks_HU.csv'),
    'Fig4C_replicates': pd.read_csv(O + 'source_data_fig4C.csv'),
    'Absent_in_rdh54_HU': G[G.absent_in_rdh54_HU.astype(bool)],
    'GO_category_tests': pd.read_csv(O + 'sga_enrichment.csv'),
    'Top200_original_method': pd.read_csv(O + 'top200_reproduced.csv'),
    'All_gene_scores': G,
}
os.makedirs('Source_Data', exist_ok=True)
with pd.ExcelWriter('Source_Data/Source_Data_Fig4_SGA.xlsx') as w:
    pd.DataFrame(README, columns=['Item', 'Description']).to_excel(w, sheet_name='README', index=False)
    for k, v in sheets.items(): v.to_excel(w, sheet_name=k, index=False)
print('Source_Data/Source_Data_Fig4_SGA.xlsx')
