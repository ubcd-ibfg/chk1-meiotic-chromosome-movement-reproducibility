"""Fig. 5F – time-lapse montage of the two representative cells (frames 1'–25').
Red frames = frames joined by a step > 1.29 µm (85th percentile of all steps) in the tracking data
(output/source_data_fig5H_steps.csv: transition t = frame t -> t+1 in 1-based minute labels).
Run after 13_cerevisiae_analysis.py:  python code/14_cerevisiae_montage.py
"""
import os as _os; _os.chdir(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))  # run from the Reproducibility folder
import logging, numpy as np, pandas as pd, tifffile, matplotlib
logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle
from scipy import ndimage
plt.rcParams.update({'font.family': ['Arial', 'Liberation Sans'], 'font.size': 7, 'pdf.fonttype': 42, 'svg.fonttype': 'none'})
PX = 0.0663; NF = 25
H = pd.read_csv('output/S_cerevisiae/source_data_fig5H_steps.csv')
ROWS = [('Control', 'C45', 'raw_data/S_cerevisiae/ctrol_3_2_C45.tif', '#2C6FB0'), ('rdh54Δ', 'C36', 'raw_data/S_cerevisiae/MAX_rdh54_02_1_C36.tif', '#E2524F')]
f, axs = plt.subplots(2, NF, figsize=(7.2, 1.15), gridspec_kw=dict(wspace=0.04, hspace=0.55))
out = []
for row, (cond, cell, tif, col) in enumerate(ROWS):
    a = tifffile.imread(tif)[:NF].astype(float)
    lo, hi = np.percentile(a, 0.5), np.percentile(a, 99.9)
    big = H[(H.cell == cell) & H.large].transition.values          # 1-based transition t: label t' -> (t+1)'
    red = sorted(set(big) | set(big + 1))                           # minute labels to frame in red
    for i in range(NF):
        ax = axs[row, i]; im = np.clip((a[i] - lo) / (hi - lo), 0, 1)
        h, w = im.shape; s = min(h, w); im = im[(h - s) // 2:(h - s) // 2 + s, (w - s) // 2:(w - s) // 2 + s]
        ax.imshow(im, cmap='gray', vmin=0, vmax=1, interpolation='nearest'); ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values(): sp.set_visible(False)
        dog = ndimage.gaussian_filter(im, 1.2) - ndimage.gaussian_filter(im, 4); y, x = np.unravel_index(dog.argmax(), dog.shape)
        ax.add_patch(Circle((x, y), 4.5, fill=False, ec=col, lw=0.6))
        lab = i + 1
        ax.text(2, 2, f"{lab}'", color='white', fontsize=5, va='top', ha='left', bbox=dict(fc='black', ec='none', pad=0.6))
        if lab in red: ax.add_patch(Rectangle((0.5, 0.5), s - 2, s - 2, fill=False, ec='#E31A1C', lw=1.1))
        out.append(dict(condition=cond, cell=cell, minute=lab, red_frame=lab in red, spot_x_px=x, spot_y_px=y))
    axs[row, 0].text(0, 1.08, ('Control  ' if row == 0 else '$rdh54\\Delta$  ') + '$TEL4L$-$GFP$', transform=axs[row, 0].transAxes, fontsize=7, va='bottom', ha='left')
    print(cond, 'red frames:', red)
# 2 µm scale bar on the first control frame
s = min(tifffile.imread(ROWS[0][2]).shape[1:]); L = 2 / PX
axs[0, 0].plot([s - 4 - L, s - 4], [s - 6, s - 6], color='white', lw=1.5, solid_capstyle='butt'); axs[0, 0].text(1.0, -0.12, '2 µm', transform=axs[0, 0].transAxes, fontsize=5, ha='right', va='top')
for e in ('pdf', 'svg'): f.savefig(f'output/S_cerevisiae/panels/Fig5F_montage.{e}', bbox_inches='tight', dpi=600)
f.savefig('output/S_cerevisiae/panels/Fig5F_montage.png', dpi=400, bbox_inches='tight')
pd.DataFrame(out).to_csv('output/S_cerevisiae/source_data_fig5F_montage.csv', index=False)
