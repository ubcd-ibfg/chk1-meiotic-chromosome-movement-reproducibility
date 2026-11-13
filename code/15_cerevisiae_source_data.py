"""Fig. 5 – step 2: assemble Source_Data_Fig5_cerevisiae.xlsx from output/*.csv (run after 13_cerevisiae_analysis.py)."""
import os as _os; _os.chdir(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))  # run from the Reproducibility folder
import os
import pandas as pd
R = ["Computed from raw_data/S_cerevisiae/control_int5_PASS_excl12_filtered.csv (46 cells) and raw_data/S_cerevisiae/mutant_int5_PASS_excl42_filtered.csv (137 cells) with code/13_cerevisiae_analysis.py.",
 "Pixel size 0.0663 µm; frame interval 1 min; 61 frames per cell. Unit of analysis: the cell.",
 "A, D: two-sided Mann–Whitney U. B: balanced accuracy, stratified 5-fold CV repeated 5 times; P from 200 label permutations (balanced accuracy). Majority-class accuracy baseline = 0.749.",
 "C: Random Forest importance (500 trees, class-balanced) for 30 per-cell features; univariate Mann–Whitney P uncorrected (no feature survives Bonferroni; minimum adjusted P = 0.34).",
 "E: % of cells with ≥1 step larger than the 85th percentile of all steps (1.29 µm); Wilson 95% CI; two-sided Fisher exact test.",
 "F: montage frames; red_frame = frame joined by a step > 1.29 µm (panel H). G, H: representative cells (control C45, field ctrol_3_2; rdh54Δ C36, field MAX_rdh54_02_1), first 25 frames, coordinates centred on the trajectory mean."]
S = {'A_D_E_per_cell': 'source_data_fig5_ADE_per_cell.csv', 'B_classifiers': 'source_data_fig5B_classifiers.csv',
     'C_feature_importance': 'source_data_fig5C_feature_importance.csv', 'E_large_displacement': 'source_data_fig5E_large_displacement.csv',
     'G_trajectories': 'source_data_fig5G_trajectories.csv', 'H_steps': 'source_data_fig5H_steps.csv', 'F_montage': 'source_data_fig5F_montage.csv'}
os.makedirs('Source_Data', exist_ok=True)
with pd.ExcelWriter('Source_Data/Source_Data_Fig5_cerevisiae.xlsx') as w:
    pd.DataFrame({'Source data – Fig. 5 (S. cerevisiae TEL4L-GFP telomere tracking)': R}).to_excel(w, sheet_name='README', index=False)
    for k, f in S.items(): pd.read_csv('output/S_cerevisiae/' + f).to_excel(w, sheet_name=k, index=False)
print('ok')
