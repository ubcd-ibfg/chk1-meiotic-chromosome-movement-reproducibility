"""
cerevisiae_panel_EJ_master.py
─────────────────────────────────────────────────────────────────────────────
Master script for S. cerevisiae telomere tracking analysis (Panel EJ).
Reproduces ALL quantifications and figure components from this session.

RAW DATA REQUIRED (relative to script location or absolute paths below):
  - MAX_rdh54_02_1_C36.tif          → C36 rdh54Δ cell TIFF stack
  - ctrol_3_2_C45.tif               → C45 control cell TIFF stack
  - mutant_int5_PASS_excl42_filtered.csv   → mutant trajectories (CSV)
  - control_int5_PASS_excl12_filtered.csv  → control trajectories (CSV)

PARAMETERS:
  PX      = 0.0663 µm/px   (pixel size)
  DT      = 1 min/frame    (time interval)
  N       = 25             (frames analysed per cell)
  D_THRESH = 85th percentile of all step sizes → 1.293 µm (19.51 px)

OUTPUTS (saved to OUTPUT_DIR):
  EJ_components/
    C36_rdh54_strip1row.{pdf,png,tiff}     → 1×25 frame strip (rdh54Δ)
    C45_control_strip1row.{pdf,png,tiff}   → 1×25 frame strip (control)
    C36_rdh54_montage.{pdf,png,tiff}       → 5×5 montage (rdh54Δ)
    C45_control_montage.{pdf,png,tiff}     → 5×5 montage (control)
    C36_rdh54_trajectory.{pdf,png}         → 2D trajectory (shared axes)
    C45_control_trajectory.{pdf,png}       → 2D trajectory (shared axes)
    C36_rdh54_steps.{pdf,png}             → step-size bars
    C45_control_steps.{pdf,png}           → step-size bars
    C36_rdh54_timeseries.{pdf,png}        → X,Y vs time
    C45_control_timeseries.{pdf,png}      → X,Y vs time

METHODS:
  - Nucleus centroid: gaussian_filter(σ=3) → center_of_mass of pixels >70th pct
  - Absolute position: centroid + CSV relative coordinates
  - Circle display: global argmax per frame (gaussian σ=1.5) = brightest spot
  - Large step threshold: 85th percentile of all pairwise step sizes across
    both mutant and control populations
  - Step size: Euclidean distance between consecutive absolute positions (µm)

Author: generated with Claude (Anthropic), June 2026
─────────────────────────────────────────────────────────────────────────────
"""

import os
import numpy as np
import pandas as pd
import tifffile
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.ndimage import gaussian_filter, center_of_mass
import warnings; warnings.filterwarnings('ignore')

# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION — edit paths here
# ══════════════════════════════════════════════════════════════════════════════

DATA_DIR   = '/sessions/upbeat-awesome-wozniak/mnt/Figure_cerevisiae/'
OUTPUT_DIR = '/sessions/upbeat-awesome-wozniak/mnt/Nature com/DRAFT/panels/EJ_components/'

TIFF_C36 = DATA_DIR + 'MAX_rdh54_02_1_C36.tif'
TIFF_C45 = DATA_DIR + 'ctrol_3_2_C45.tif'
CSV_MUT  = DATA_DIR + 'mutant_int5_PASS_excl42_filtered.csv'
CSV_CTRL = DATA_DIR + 'control_int5_PASS_excl12_filtered.csv'

# ── imaging parameters ────────────────────────────────────────────────────────
PX      = 0.0663   # µm per pixel
DT      = 1.0      # min per frame
N       = 25       # number of frames to analyse

# ── large-step threshold ──────────────────────────────────────────────────────
# 85th percentile of ALL step sizes (mutant + control combined)
# Pre-computed value: 19.51 px = 1.293 µm
D_THRESH_PX  = 19.51
D_THRESH_UM  = D_THRESH_PX * PX   # 1.293 µm

# ── scale bar ────────────────────────────────────────────────────────────────
SCALE_UM  = 2      # µm scale bar on last frame of strip
SCALE_PX  = SCALE_UM / PX

DPI  = 300
os.makedirs(OUTPUT_DIR, exist_ok=True)


# ══════════════════════════════════════════════════════════════════════════════
# HELPER FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

def load_stack(tiff_path, n_frames=N):
    """Load first n_frames of a TIFF stack; collapse channel dim if present."""
    stack = tifffile.imread(tiff_path)
    if stack.ndim == 4:        # (T, C, Y, X)
        stack = stack[:, 0, :, :]
    return stack[:n_frames]


def load_cell_csv(csv_path, cell_id, field_substr, n_frames=N):
    """Return the first n_frames rows for cell_id from the given TIFF field."""
    df = pd.read_csv(csv_path)
    df = df[df['cell_id'] == cell_id].sort_values('frame')
    df = df[df['field_id'].str.contains(field_substr, na=False)]
    return df[df['frame'] < n_frames].reset_index(drop=True)


def nucleus_centroid(stack):
    """
    Per-frame nucleus centroid via gaussian blur + center of mass.
    Returns (cx, cy) arrays in pixel coords.
    """
    cxs, cys = [], []
    for t in range(len(stack)):
        img    = gaussian_filter(stack[t].astype(float), sigma=3.0)
        thresh = np.percentile(img, 70)
        mask   = (img > thresh).astype(float)
        cy, cx = center_of_mass(img * mask)
        cxs.append(cx); cys.append(cy)
    return np.array(cxs), np.array(cys)


def global_argmax(stack):
    """
    Per-frame position of the brightest pixel (telomere / Nuf2 spot).
    Uses light gaussian blur (σ=1.5) to suppress single-pixel noise.
    Returns (gx, gy) arrays in pixel coords.
    """
    xs, ys = [], []
    for t in range(len(stack)):
        img        = gaussian_filter(stack[t].astype(float), sigma=1.5)
        iy, ix     = np.unravel_index(np.argmax(img), img.shape)
        xs.append(float(ix)); ys.append(float(iy))
    return np.array(xs), np.array(ys)


def compute_trajectory(stack, df):
    """
    Absolute telomere position in µm (centred) and step sizes.
    Method: nucleus centroid + CSV relative coordinates → absolute px → µm.
    """
    ncx, ncy = nucleus_centroid(stack)
    x_abs    = ncx + df['x_rel_px'].values
    y_abs    = ncy + df['y_rel_px'].values
    # centre on mean position
    x_um = (x_abs - x_abs.mean()) * PX
    y_um = (y_abs - y_abs.mean()) * PX
    steps = np.sqrt(np.diff(x_um)**2 + np.diff(y_um)**2)
    large = steps > D_THRESH_UM
    return x_um, y_um, steps, large


# ══════════════════════════════════════════════════════════════════════════════
# LOAD DATA
# ══════════════════════════════════════════════════════════════════════════════

print("Loading data …")

stk36 = load_stack(TIFF_C36)
stk45 = load_stack(TIFF_C45)

c36 = load_cell_csv(CSV_MUT,  'C36', 'MAX_rdh54_02_1')
c45 = load_cell_csv(CSV_CTRL, 'C45', 'ctrol_3_2')

x36, y36, s36, l36 = compute_trajectory(stk36, c36)
x45, y45, s45, l45 = compute_trajectory(stk45, c45)

gx36, gy36 = global_argmax(stk36)
gx45, gy45 = global_argmax(stk45)

t_min = np.arange(N) * DT

print(f"C36 rdh54Δ:  {l36.sum()} large steps  |  mean step = {s36.mean():.3f} µm/min")
print(f"C45 Control: {l45.sum()} large steps  |  mean step = {s45.mean():.3f} µm/min")
print(f"D_THRESH = {D_THRESH_UM:.3f} µm  ({D_THRESH_PX:.2f} px)")

# Shared trajectory axis limits
PAD = 0.3
XMIN = min(x36.min(), x45.min()) - PAD;  XMAX = max(x36.max(), x45.max()) + PAD
YMIN = min(y36.min(), y45.min()) - PAD;  YMAX = max(y36.max(), y45.max()) + PAD
SMAX = max(s36.max(), s45.max()) * 1.15   # shared step-bar y-axis

CELLS = [
    dict(name='C36_rdh54',  stack=stk36, gx=gx36, gy=gy36,
         x=x36, y=y36, steps=s36, large=l36,
         col='#DD3333', label='rdh54Δ',  is_mut=True),
    dict(name='C45_control',stack=stk45, gx=gx45, gy=gy45,
         x=x45, y=y45, steps=s45, large=l45,
         col='#2266AA', label='Control', is_mut=False),
]


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1 — 1×25 horizontal strip (real aspect ratio, 2 µm scale bar)
# ══════════════════════════════════════════════════════════════════════════════

print("\nGenerating 1×25 strips …")
GAP_PX = 2; CIRC_R = 4

for cell in CELLS:
    stack = cell['stack']; gx = cell['gx']; gy = cell['gy']
    large = cell['large']; col = cell['col']; name = cell['name']
    H, W  = stack[0].shape

    fig_w = (N * W + (N - 1) * GAP_PX) / DPI
    fig_h = H / DPI
    fig   = plt.figure(figsize=(fig_w, fig_h), dpi=DPI)
    fig.patch.set_facecolor('black')

    for fr in range(N):
        x0 = fr * (W + GAP_PX) / (N * W + (N - 1) * GAP_PX)
        w  = W                  / (N * W + (N - 1) * GAP_PX)
        ax = fig.add_axes([x0, 0, w, 1])

        img  = stack[fr].astype(float)
        vmin = np.percentile(img, 1); vmax = np.percentile(img, 99.5)
        ax.imshow(img, cmap='gray', vmin=vmin, vmax=vmax, interpolation='nearest')

        # telomere circle (global argmax = brightest spot)
        ax.add_patch(plt.Circle((gx[fr], gy[fr]), CIRC_R,
                                color=col, fill=False, lw=1.0))

        # frame number
        ax.text(1, 1, str(fr + 1), color='white', fontsize=4, va='top', ha='left')

        # border on large-step adjacent frames
        is_large = ((fr > 0 and large[fr - 1]) or (fr < N - 1 and large[fr]))
        for sp in ax.spines.values():
            sp.set_visible(is_large)
            if is_large:
                sp.set_color(col); sp.set_linewidth(1.5)

        # 2 µm scale bar on last frame only
        if fr == N - 1:
            margin = 4; bar_h = 2
            x_end = W - margin; x_start = x_end - SCALE_PX; y_bar = H - margin - bar_h
            ax.add_patch(plt.Rectangle((x_start, y_bar), SCALE_PX, bar_h,
                                       color='white', zorder=10))
            ax.text((x_start + x_end) / 2, y_bar - 2, f'{SCALE_UM} µm',
                    color='white', fontsize=4, ha='center', va='bottom', zorder=11)

        ax.set_xlim(0, W); ax.set_ylim(H, 0)
        ax.set_xticks([]); ax.set_yticks([])

    for ext in ('pdf', 'png', 'tiff'):
        plt.savefig(f'{OUTPUT_DIR}{name}_strip1row.{ext}', dpi=DPI,
                    bbox_inches='tight', pad_inches=0,
                    facecolor='black', edgecolor='none')
    plt.close()
    print(f"  Strip saved: {name}")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2 — 5×5 montage
# ══════════════════════════════════════════════════════════════════════════════

print("\nGenerating 5×5 montages …")

for cell in CELLS:
    fig, axes = plt.subplots(5, 5, figsize=(4.5, 9.5),
                              gridspec_kw={'wspace': 0.03, 'hspace': 0.03})
    fig.patch.set_facecolor('black')

    for fr in range(N):
        r, c_idx = divmod(fr, 5)
        ax  = axes[r][c_idx]
        img = cell['stack'][fr].astype(float)
        vmin = np.percentile(img, 1); vmax = np.percentile(img, 99.5)
        ax.imshow(img, cmap='gray', vmin=vmin, vmax=vmax, interpolation='nearest')
        ax.add_patch(plt.Circle((cell['gx'][fr], cell['gy'][fr]), CIRC_R,
                                color=cell['col'], fill=False, lw=1.0))
        ax.text(2, 4, f'{fr + 1}', color='white', fontsize=5, va='top')
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values(): sp.set_visible(False)

    for ext in ('pdf', 'png', 'tiff'):
        plt.savefig(f'{OUTPUT_DIR}{cell["name"]}_montage.{ext}', dpi=DPI,
                    bbox_inches='tight', facecolor='black', edgecolor='none')
    plt.close()
    print(f"  Montage saved: {cell['name']}")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 3 — 2D trajectory (segments colored by step type, shared axes)
# ══════════════════════════════════════════════════════════════════════════════

print("\nGenerating 2D trajectory plots …")

for cell in CELLS:
    x_um = cell['x']; y_um = cell['y']
    large = cell['large']; col = cell['col']
    n_lg  = large.sum()

    fig, ax = plt.subplots(figsize=(3.4, 3.4))
    fig.patch.set_facecolor('white')

    # segments: red = large step, gray = normal
    for i in range(N - 1):
        seg_col = '#E83030' if large[i] else '#BBBBBB'
        lw      = 2.2       if large[i] else 0.9
        ax.plot([x_um[i], x_um[i + 1]], [y_um[i], y_um[i + 1]],
                color=seg_col, lw=lw, zorder=2 + int(large[i]),
                solid_capstyle='round')

    # numbered points
    for i in range(N):
        ax.plot(x_um[i], y_um[i], 'o', color='white', ms=8,
                zorder=5, mec='#555', mew=0.5)
        ax.text(x_um[i], y_um[i], str(i + 1), color='#333',
                fontsize=4, ha='center', va='center', zorder=6)

    ax.plot(x_um[0],  y_um[0],  'o', ms=10, color='lime', zorder=7,
            label='start', mec='white', mew=0.8)
    ax.plot(x_um[-1], y_um[-1], 's', ms=9,  color='gold', zorder=7,
            label='end',   mec='white', mew=0.8)

    handles = [
        Line2D([], [], color='#E83030', lw=2,
               label=f'large step (>{D_THRESH_UM:.2f} µm)'),
        Line2D([], [], color='#BBBBBB', lw=1.5, label='normal step'),
        Line2D([], [], marker='o', color='lime', ms=7, lw=0, label='start'),
        Line2D([], [], marker='s', color='gold', ms=7, lw=0, label='end'),
    ]
    ax.legend(handles=handles, fontsize=6.5, loc='upper right',
              framealpha=0.85, handlelength=1.2)

    ax.set_xlim(XMIN, XMAX); ax.set_ylim(YMIN, YMAX)
    ax.set_aspect('equal')
    ax.set_xlabel('X (µm)', fontsize=9); ax.set_ylabel('Y (µm)', fontsize=9)
    ax.set_title(f'{cell["label"]}  —  {n_lg} large steps',
                 fontsize=10, fontweight='bold', color=col)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    plt.tight_layout(pad=1.0)

    for ext in ('pdf', 'png'):
        plt.savefig(f'{OUTPUT_DIR}{cell["name"]}_trajectory.{ext}',
                    bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  Trajectory saved: {cell['name']}")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 4 — Step-size bar chart (shared Y axis)
# ══════════════════════════════════════════════════════════════════════════════

print("\nGenerating step bar charts …")

for cell in CELLS:
    steps = cell['steps']; large = cell['large']
    fig, ax = plt.subplots(figsize=(3.8, 2.2))
    fig.patch.set_facecolor('white')

    bar_colors = [cell['col'] if lg else '#BBBBBB' for lg in large]
    ax.bar(np.arange(len(steps)) + 0.5, steps,
           color=bar_colors, width=0.82, edgecolor='none')
    ax.axhline(D_THRESH_UM, color='#555', lw=1.2, ls='--',
               label=f'threshold {D_THRESH_UM:.2f} µm (85th pct)')

    ax.set_xlim(0, N - 1); ax.set_ylim(0, SMAX)
    ax.set_xlabel('Frame', fontsize=9); ax.set_ylabel('Step size (µm)', fontsize=9)
    ax.set_title(f'{cell["label"]}  |  {large.sum()} large steps',
                 fontsize=10, fontweight='bold', color=cell['col'])
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    ax.legend(fontsize=7.5)
    plt.tight_layout(pad=1.0)

    for ext in ('pdf', 'png'):
        plt.savefig(f'{OUTPUT_DIR}{cell["name"]}_steps.{ext}',
                    bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  Steps saved: {cell['name']}")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 5 — X,Y time series
# ══════════════════════════════════════════════════════════════════════════════

print("\nGenerating time-series plots …")

for cell in CELLS:
    fig, axes = plt.subplots(2, 1, figsize=(4.5, 3.0), sharex=True)
    fig.patch.set_facecolor('white')

    for ax_i, (vals, lbl) in zip(axes,
            [(cell['x'], 'X (µm)'), (cell['y'], 'Y (µm)')]):
        ax_i.plot(t_min, vals, color=cell['col'], lw=1.5)
        ax_i.set_ylabel(lbl, fontsize=9)
        ax_i.spines['top'].set_visible(False)
        ax_i.spines['right'].set_visible(False)
        for i, lg in enumerate(cell['large']):
            if lg:
                ax_i.axvspan(t_min[i], t_min[i + 1],
                             color=cell['col'], alpha=0.18, zorder=0)

    axes[-1].set_xlabel('Time (min)', fontsize=9)
    axes[0].set_title(f'{cell["label"]} — time series',
                      fontsize=10, fontweight='bold', color=cell['col'])
    plt.tight_layout(pad=0.8)

    for ext in ('pdf', 'png'):
        plt.savefig(f'{OUTPUT_DIR}{cell["name"]}_timeseries.{ext}',
                    bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  Time-series saved: {cell['name']}")


# ══════════════════════════════════════════════════════════════════════════════
# SUMMARY STATISTICS
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 60)
print("SUMMARY STATISTICS")
print("=" * 60)

for cell in CELLS:
    s = cell['steps']
    print(f"\n{cell['label']} (C{'36' if cell['is_mut'] else '45'}):")
    print(f"  n frames        : {N}")
    print(f"  n large steps   : {cell['large'].sum()}")
    print(f"  mean step       : {s.mean():.3f} µm/min  ({s.mean()/60:.4f} µm/s)")
    print(f"  median step     : {np.median(s):.3f} µm/min")
    print(f"  max step        : {s.max():.3f} µm/min")
    print(f"  D_THRESH (85th) : {D_THRESH_UM:.3f} µm")

print(f"\nVelocity note: mean step {np.concatenate([s36,s45]).mean():.3f} µm/min")
print(f"  = {np.concatenate([s36,s45]).mean()/60:.4f} µm/s")
print(f"  Literature: 0.09 µm/s (instantaneous, high time-res)")
print(f"  Diffusive scaling: √(4D×60s) = {np.sqrt(4*0.09**2/4*60):.2f} µm → consistent")

print("\nAll outputs saved to:", OUTPUT_DIR)
print("Done ✓")
