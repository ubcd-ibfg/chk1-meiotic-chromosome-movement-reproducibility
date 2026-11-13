"""Fig. 1F z-score heat map as a pure vector panel (cells drawn as vector quads, text ≥ 6 pt). Reads output/fig1f_zscores.csv (01)."""
from common import *
import logging; logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)
plt.rcParams.update({'font.family': ['Arial', 'Liberation Sans'], 'font.size': 7})
Z = pd.read_csv(os.path.join(OUT, 'fig1f_zscores.csv'), index_col=0)
G = ORDER[1:]; Z = Z.loc[G]; LB = ['Area', 'Major axis', 'Minor axis', 'Circularity', 'Convexity', 'Velocity']
f, ax = plt.subplots(figsize=(2.95, 2.7))
pc = ax.pcolormesh(np.arange(7) - 0.5, np.arange(14) - 0.5, Z.values, cmap='RdBu_r', vmin=-4, vmax=4, edgecolors='none')
for i in range(13):
    for j in range(6):
        v = Z.values[i, j]; ax.text(j, i, f'{v:.2f}', ha='center', va='center', fontsize=6, color='white' if abs(v) > 2.5 else 'black')
ax.set_xlim(-0.5, 5.5); ax.set_ylim(12.5, -0.5)
ax.set_xticks(range(6)); ax.set_xticklabels(LB, rotation=40, ha='right'); ax.set_yticks(range(13)); ax.set_yticklabels(G, fontstyle='italic')
for s in ax.spines.values(): s.set_visible(False)
for g0, g1, c in [(0, 0, '#E8743B'), (1, 4, '#D6338A'), (5, 10, '#6A5ACD'), (11, 11, '#2A9D6F'), (12, 12, '#D9A21B')]:
    ax.plot([-2.75, -2.75], [g0 - 0.45, g1 + 0.45], color=c, lw=4, clip_on=False, solid_capstyle='butt')
cb = f.colorbar(pc, ax=ax, shrink=0.75, pad=0.03, ticks=range(-4, 5)); cb.set_label('z-score vs control'); cb.outline.set_linewidth(0.5); cb.solids.set_rasterized(False)
ax.set_title('Phenotypic deviation from control', fontsize=7)
os.makedirs(os.path.join(OUT, 'vector_panels'), exist_ok=True)
for e in ('pdf', 'svg'): f.savefig(os.path.join(OUT, 'vector_panels', f'Fig1F_zscores.{e}'), bbox_inches='tight')
print('ok')
