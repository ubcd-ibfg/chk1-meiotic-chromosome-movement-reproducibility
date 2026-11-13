#!/usr/bin/env python3
"""
S. cerevisiae meiotic telomere trajectory analysis
Control vs rdh54Δ

Two complementary approaches:
  1. Machine learning classification (Random Forest + 8 classifiers)
     — Feature extraction from per-cell trajectory descriptors
     — Stratified 5-fold cross-validation
     — Feature importance ranking

  2. Large-displacement event analysis (Motif C-equivalent search)
     — Define "large steps" as displacements > 85th percentile
     — Characterise "return events": large step followed by return
       toward origin within 2-4 frames (analogous to Motif C reset)
     — Compare frequencies between control and rdh54Δ
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis, QuadraticDiscriminantAnalysis
from sklearn.svm import SVC
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import confusion_matrix
import warnings
warnings.filterwarnings('ignore')

DT = 5.0   # min per frame

# ─────────────────────────────────────────────────────────────────────────────
# 1. LOAD DATA
# ─────────────────────────────────────────────────────────────────────────────
ctrl = pd.read_csv('/sessions/upbeat-awesome-wozniak/mnt/uploads/control_int5_PASS_excl12_filtered.csv')
mut  = pd.read_csv('/sessions/upbeat-awesome-wozniak/mnt/uploads/mutant_int5_PASS_excl42_filtered.csv')
ctrl['uid']   = ctrl['field_id'] + '__' + ctrl['cell_id']
mut['uid']    = mut['field_id']  + '__' + mut['cell_id']
ctrl['label'] = 0   # 0 = control
mut['label']  = 1   # 1 = rdh54Δ
ctrl = ctrl.sort_values(['uid','frame']).reset_index(drop=True)
mut  = mut.sort_values(['uid','frame']).reset_index(drop=True)
df_all = pd.concat([ctrl, mut], ignore_index=True)

# ─────────────────────────────────────────────────────────────────────────────
# 2. FEATURE EXTRACTION (per-cell, trajectory-level)
# ─────────────────────────────────────────────────────────────────────────────

# Compute global step threshold for large-displacement detection
all_steps = []
for uid, sub in df_all.groupby('uid'):
    sub = sub.sort_values('frame')
    x, y = sub['x_rel_px'].values, sub['y_rel_px'].values
    s = np.sqrt(np.diff(x)**2 + np.diff(y)**2)
    all_steps.extend(s[s > 0].tolist())
all_steps = np.array(all_steps)
THRESH_85 = np.percentile(all_steps, 85)
THRESH_90 = np.percentile(all_steps, 90)
THRESH_95 = np.percentile(all_steps, 95)
print(f"Global step thresholds: 85th={THRESH_85:.1f}px, 90th={THRESH_90:.1f}px, 95th={THRESH_95:.1f}px")


def extract_features(df):
    records = []
    for uid, sub in df.groupby('uid'):
        sub  = sub.sort_values('frame')
        lbl  = sub['label'].iloc[0]
        x, y = sub['x_rel_px'].values, sub['y_rel_px'].values
        n    = len(x)
        dx, dy = np.diff(x), np.diff(y)
        step   = np.sqrt(dx**2 + dy**2)
        sp     = sub['speed_px_frame'].dropna().values

        # ── Basic kinematics ────────────────────────────────────────────
        mean_step  = step.mean()
        std_step   = step.std()
        max_step   = step.max()
        p25_step   = np.percentile(step, 25)
        p75_step   = np.percentile(step, 75)
        iqr_step   = p75_step - p25_step

        # ── Displacement distribution shape ─────────────────────────────
        # Kurtosis and skewness of step distribution
        kurt_step  = stats.kurtosis(step)
        skew_step  = stats.skew(step)
        # Non-Gaussian parameter (van Hove)
        r2 = dx**2 + dy**2
        ng = np.mean(r2**2) / (2 * np.mean(r2)**2) - 1 if np.mean(r2) > 0 else 0

        # ── Directional statistics ──────────────────────────────────────
        angles = np.arctan2(dy, dx)
        # VACF at lags 1, 2, 3
        v2 = np.mean(dx**2 + dy**2)
        if v2 > 0:
            vacf1 = np.mean(dx[:-1]*dx[1:] + dy[:-1]*dy[1:]) / v2
            vacf2 = np.mean(dx[:-2]*dx[2:] + dy[:-2]*dy[2:]) / v2 if len(dx)>2 else 0
            vacf3 = np.mean(dx[:-3]*dx[3:] + dy[:-3]*dy[3:]) / v2 if len(dx)>3 else 0
        else:
            vacf1 = vacf2 = vacf3 = 0

        # Directional reversal rate (fraction of steps reversing direction)
        if len(dx) > 1:
            dots = dx[:-1]*dx[1:] + dy[:-1]*dy[1:]
            rev_rate = np.mean(dots < 0)
        else:
            rev_rate = 0

        # ── MSD and anomalous diffusion ─────────────────────────────────
        lags = np.arange(1, min(11, n))
        msds = [np.mean((x[l:]-x[:-l])**2 + (y[l:]-y[:-l])**2) for l in lags]
        # Anomalous diffusion exponent α
        if len(lags) >= 3 and all(m > 0 for m in msds):
            alpha = np.polyfit(np.log(lags * DT), np.log(msds), 1)[0]
        else:
            alpha = 1.0
        # MSD plateau (ratio MSD_10/MSD_1 — low = confined)
        msd_ratio = msds[-1] / msds[0] if len(msds) >= 2 and msds[0] > 0 else 1.0
        msd_lag1  = msds[0]
        msd_lag5  = msds[4] if len(msds) > 4 else np.nan

        # ── Spatial statistics ──────────────────────────────────────────
        rg   = np.sqrt(np.var(x) + np.var(y))    # radius of gyration
        span = np.sqrt((x.max()-x.min())**2 + (y.max()-y.min())**2)
        # Confinement: net displacement / total path
        net_disp = np.sqrt((x[-1]-x[0])**2 + (y[-1]-y[0])**2)
        cum_dist = step.sum()
        conf     = net_disp / cum_dist if cum_dist > 0 else 0

        # ── Large displacement events ────────────────────────────────────
        large85   = np.sum(step >= THRESH_85)
        large90   = np.sum(step >= THRESH_90)
        large95   = np.sum(step >= THRESH_95)
        frac_l85  = large85 / len(step)
        frac_l95  = large95 / len(step)

        # ── Return events (Motif C-equivalent) ──────────────────────────
        # Definition: a "large step" (>85th pct) followed within W frames
        # by a return that brings the telomere back within RETURN_FRAC
        # of its pre-large-step position
        return2, return3, return4 = 0, 0, 0
        for i in range(len(step)):
            if step[i] >= THRESH_85:
                x_pre, y_pre = x[i], y[i]
                # position after jump is x[i+1], y[i+1]
                for W, slot in [(2, 0), (3, 1), (4, 2)]:
                    if i + 1 + W < n:   # strict < to avoid out-of-bounds
                        x_check = x[i + 1 + W]; y_check = y[i + 1 + W]
                        dist_to_pre  = np.sqrt((x_check - x_pre)**2 + (y_check - y_pre)**2)
                        dist_of_jump = step[i]
                        if dist_to_pre < 0.5 * dist_of_jump:
                            if slot == 0: return2 += 1
                            elif slot == 1: return3 += 1
                            elif slot == 2: return4 += 1

        # Also: large step followed by large step in OPPOSITE direction (direct reversal)
        reversals_large = 0
        for i in range(len(step)-1):
            if step[i] >= THRESH_85 and step[i+1] >= THRESH_85:
                dot = dx[i]*dx[i+1] + dy[i]*dy[i+1]
                if dot < 0:
                    reversals_large += 1

        # ── Temporal autocorrelation of step magnitude ───────────────────
        # Does the step size itself have structure? (periodic events?)
        acf_step = np.correlate(step - step.mean(), step - step.mean(), mode='full')
        acf_step = acf_step[len(step)-1:] / acf_step[len(step)-1]   # normalise
        acf_step_lag1 = acf_step[1] if len(acf_step) > 1 else 0
        acf_step_lag3 = acf_step[3] if len(acf_step) > 3 else 0

        records.append({
            'uid': uid, 'label': lbl, 'n_frames': n,
            'mean_step': mean_step, 'std_step': std_step, 'max_step': max_step,
            'p25_step': p25_step, 'p75_step': p75_step, 'iqr_step': iqr_step,
            'kurt_step': kurt_step, 'skew_step': skew_step, 'ng_param': ng,
            'vacf1': vacf1, 'vacf2': vacf2, 'vacf3': vacf3,
            'rev_rate': rev_rate,
            'alpha': alpha, 'msd_ratio': msd_ratio, 'msd_lag1': msd_lag1,
            'rg': rg, 'span': span, 'conf': conf,
            'large85': large85, 'large90': large90, 'large95': large95,
            'frac_l85': frac_l85, 'frac_l95': frac_l95,
            'return_2fr': return2, 'return_3fr': return3, 'return_4fr': return4,
            'rev_large': reversals_large,
            'acf_step_lag1': acf_step_lag1, 'acf_step_lag3': acf_step_lag3,
        })
    return pd.DataFrame(records)


print("Extracting features…")
feat_ctrl = extract_features(ctrl)
feat_mut  = extract_features(mut)
feat_all  = pd.concat([feat_ctrl, feat_mut], ignore_index=True)
print(f"  Control: {len(feat_ctrl)} cells | rdh54Δ: {len(feat_mut)} cells")
print(f"  Features: {len(feat_all.columns)-3}")


# ─────────────────────────────────────────────────────────────────────────────
# 3. ML CLASSIFICATION
# ─────────────────────────────────────────────────────────────────────────────
FEATURE_COLS = [c for c in feat_all.columns if c not in ('uid','label','n_frames')]
X = feat_all[FEATURE_COLS].fillna(0).values
y = feat_all['label'].values

classifiers = {
    'Logistic Reg.'  : Pipeline([('sc', StandardScaler()), ('clf', LogisticRegression(max_iter=500))]),
    'Decision Tree'  : DecisionTreeClassifier(random_state=42),
    'Naive Bayes'    : GaussianNB(),
    'KNN'            : Pipeline([('sc', StandardScaler()), ('clf', KNeighborsClassifier(n_neighbors=5))]),
    'SVM'            : Pipeline([('sc', StandardScaler()), ('clf', SVC(kernel='rbf'))]),
    'Random Forest'  : RandomForestClassifier(n_estimators=200, random_state=42),
    'LDA'            : Pipeline([('sc', StandardScaler()), ('clf', LinearDiscriminantAnalysis())]),
    'QDA'            : Pipeline([('sc', StandardScaler()), ('clf', QuadraticDiscriminantAnalysis())]),
}

skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
results = {}
print("\nClassifier benchmark (5-fold CV):")
for name, clf in classifiers.items():
    scores = cross_val_score(clf, X, y, cv=skf, scoring='accuracy')
    results[name] = scores
    print(f"  {name:<18}: {scores.mean():.3f} ± {scores.std():.3f}   (min={scores.min():.3f}, max={scores.max():.3f})")

# Best classifier: Random Forest → feature importances
rf = RandomForestClassifier(n_estimators=500, random_state=42)
rf.fit(X, y)
importances = pd.Series(rf.feature_importances_, index=FEATURE_COLS).sort_values(ascending=False)
print("\nTop 15 feature importances (Random Forest):")
for feat, imp in importances.head(15).items():
    print(f"  {feat:<22}: {imp:.4f}")


# ─────────────────────────────────────────────────────────────────────────────
# 4. LARGE DISPLACEMENT EVENT ANALYSIS
# ─────────────────────────────────────────────────────────────────────────────
print("\n── Large displacement event analysis ──────────────────────────────────")
for col in ['large85','large90','large95','frac_l85','return_2fr','return_3fr','return_4fr','rev_large','kurt_step']:
    c = feat_ctrl[col].values; m = feat_mut[col].values
    _, pv = stats.mannwhitneyu(c, m, alternative='two-sided')
    sig = '***' if pv<0.001 else '**' if pv<0.01 else '*' if pv<0.05 else 'ns'
    fold = m.mean()/c.mean() if c.mean()>0 else np.nan
    print(f"  {col:<20}: ctrl={c.mean():.3f}±{(c.std()/np.sqrt(len(c))):.3f}, "
          f"mut={m.mean():.3f}±{(m.std()/np.sqrt(len(m))):.3f}, p={pv:.4f} {sig}, "
          f"fold={fold:.2f}x")

# ── Return event characterisation ────────────────────────────────────────────
print("\n── Return event characterisation (% of large steps that return) ──────")
for W, col in [(2,'return_2fr'), (3,'return_3fr'), (4,'return_4fr')]:
    c_ret  = feat_ctrl[col].values
    c_lge  = feat_ctrl['large85'].values
    m_ret  = feat_mut[col].values
    m_lge  = feat_mut['large85'].values
    c_pct  = np.divide(c_ret, c_lge, out=np.zeros_like(c_ret, dtype=float), where=c_lge>0)
    m_pct  = np.divide(m_ret, m_lge, out=np.zeros_like(m_ret, dtype=float), where=m_lge>0)
    _, pv  = stats.mannwhitneyu(c_pct, m_pct, alternative='two-sided')
    sig = '***' if pv<0.001 else '**' if pv<0.01 else '*' if pv<0.05 else 'ns'
    print(f"  W={W} frames ({W*DT:.0f} min): ctrl={c_pct.mean()*100:.1f}%, "
          f"mut={m_pct.mean()*100:.1f}%, p={pv:.4f} {sig}")


# ─────────────────────────────────────────────────────────────────────────────
# 5. FIGURES
# ─────────────────────────────────────────────────────────────────────────────
CCOL = '#444444'; MCOL = '#CC2222'
rng  = np.random.default_rng(7)

fig = plt.figure(figsize=(20, 16))
fig.patch.set_facecolor('white')
gs  = gridspec.GridSpec(4, 4, figure=fig, hspace=0.60, wspace=0.40)

# ── Panel A: Classifier benchmark ────────────────────────────────────────────
ax = fig.add_subplot(gs[0, :2])
ax.set_facecolor('white')
names = list(results.keys())
means = [results[n].mean() for n in names]
sems  = [results[n].std() / np.sqrt(5) for n in names]
colors = ['#CC2222' if n == 'Random Forest' else '#AAAAAA' for n in names]
bars = ax.barh(names, means, xerr=sems, color=colors, edgecolor='#333333',
               linewidth=0.8, capsize=3, height=0.6)
ax.axvline(0.5, color='#333333', lw=1, linestyle='--', label='chance (50%)')
ax.set_xlabel('Accuracy (5-fold CV)', fontsize=8)
ax.set_title('A   ML classifier benchmark (Control vs rdh54Δ)', fontweight='bold', loc='left', fontsize=9)
ax.tick_params(labelsize=7.5); ax.set_xlim(0.3, 0.9)
ax.legend(fontsize=7); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
# Annotate best
best_idx = np.argmax(means)
ax.text(means[best_idx]+sems[best_idx]+0.01, best_idx, f'{means[best_idx]:.3f}', va='center', fontsize=7.5, color='#CC2222', fontweight='bold')

# ── Panel B: Feature importances ─────────────────────────────────────────────
ax = fig.add_subplot(gs[0, 2:])
ax.set_facecolor('white')
top15 = importances.head(15)
fcols = ['#CC4400' if 'return' in f or 'large' in f or 'rev_l' in f or 'kurt' in f
         else '#6688AA' for f in top15.index]
ax.barh(range(15), top15.values[::-1], color=fcols[::-1], edgecolor='#333333', linewidth=0.7, height=0.7)
ax.set_yticks(range(15))
ax.set_yticklabels(top15.index[::-1], fontsize=7.5)
ax.set_xlabel('Feature importance (RF, n=500 trees)', fontsize=8)
ax.set_title('B   Random Forest feature importances\n(orange = large-displacement features)', fontweight='bold', loc='left', fontsize=9)
ax.tick_params(labelsize=7.5); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

# ── Panel C: Large steps per cell ────────────────────────────────────────────
ax = fig.add_subplot(gs[1, :2])
ax.set_facecolor('white')
data_lge = [feat_ctrl['large85'].values, feat_mut['large85'].values]
bp = ax.boxplot(data_lge, patch_artist=True, widths=0.5,
                medianprops=dict(color='black',lw=2),
                flierprops=dict(marker='o', ms=3, alpha=0.4))
for p, c in zip(bp['boxes'], [CCOL, MCOL]): p.set_facecolor(c); p.set_alpha(0.65)
for xi, (d, c) in enumerate(zip(data_lge, [CCOL, MCOL]), 1):
    ax.scatter(np.full(len(d),xi)+rng.uniform(-0.2,0.2,len(d)), d, s=12, alpha=0.5, color=c, zorder=3)
_, pv = stats.mannwhitneyu(data_lge[0], data_lge[1], alternative='two-sided')
sig = '***' if pv<0.001 else '**' if pv<0.01 else '*' if pv<0.05 else 'ns'
ax.text(1.5, max(data_lge[0].max(), data_lge[1].max())+0.3, sig+f' (p={pv:.3f})', ha='center', fontsize=10)
ax.set_xticks([1,2]); ax.set_xticklabels(['Control', r'$rdh54\Delta$'], fontsize=9)
ax.set_ylabel('Large steps per cell\n(> 85th percentile)', fontsize=8)
ax.set_title('C   Large displacement events per cell', fontweight='bold', loc='left', fontsize=9)
ax.tick_params(labelsize=7); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

# ── Panel D: Kurtosis of step distribution ────────────────────────────────────
ax = fig.add_subplot(gs[1, 2:])
ax.set_facecolor('white')
data_kurt = [feat_ctrl['kurt_step'].values, feat_mut['kurt_step'].values]
bp = ax.boxplot(data_kurt, patch_artist=True, widths=0.5,
                medianprops=dict(color='black',lw=2),
                flierprops=dict(marker='o', ms=3, alpha=0.4))
for p, c in zip(bp['boxes'], [CCOL, MCOL]): p.set_facecolor(c); p.set_alpha(0.65)
for xi, (d, c) in enumerate(zip(data_kurt, [CCOL, MCOL]), 1):
    ax.scatter(np.full(len(d),xi)+rng.uniform(-0.2,0.2,len(d)), d, s=12, alpha=0.5, color=c, zorder=3)
_, pv_k = stats.mannwhitneyu(data_kurt[0], data_kurt[1], alternative='two-sided')
sig_k = '***' if pv_k<0.001 else '**' if pv_k<0.01 else '*' if pv_k<0.05 else 'ns'
ax.text(1.5, max(d.max() for d in data_kurt)+0.3, sig_k+f' (p={pv_k:.3f})', ha='center', fontsize=10)
ax.axhline(0, color='#999', lw=0.8, linestyle='--', label='Gaussian')
ax.set_xticks([1,2]); ax.set_xticklabels(['Control', r'$rdh54\Delta$'], fontsize=9)
ax.set_ylabel('Kurtosis of step distribution\n(excess; 0 = Gaussian)', fontsize=8)
ax.set_title('D   Step distribution kurtosis\n(heavy tail = rare large events)', fontweight='bold', loc='left', fontsize=9)
ax.legend(fontsize=7); ax.tick_params(labelsize=7); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

# ── Panel E: Return events (Motif C-equivalent) ───────────────────────────────
ax = fig.add_subplot(gs[2, :2])
ax.set_facecolor('white')
windows = [2, 3, 4]
labels_w = ['10 min', '15 min', '20 min']
x_pos = np.arange(len(windows))
w_ctrl, w_mut, pvals = [], [], []
for W, col in zip(windows, ['return_2fr','return_3fr','return_4fr']):
    c_ret = feat_ctrl[col].values; c_lge = feat_ctrl['large85'].values
    m_ret = feat_mut[col].values;  m_lge = feat_mut['large85'].values
    c_pct = np.divide(c_ret, c_lge, out=np.zeros_like(c_ret,dtype=float), where=c_lge>0)
    m_pct = np.divide(m_ret, m_lge, out=np.zeros_like(m_ret,dtype=float), where=m_lge>0)
    w_ctrl.append(c_pct.mean()*100); w_mut.append(m_pct.mean()*100)
    _, pv = stats.mannwhitneyu(c_pct, m_pct, alternative='two-sided')
    pvals.append(pv)
width = 0.35
ax.bar(x_pos-width/2, w_ctrl, width, color=CCOL, alpha=0.75, edgecolor='#333', linewidth=0.8, label='Control')
ax.bar(x_pos+width/2, w_mut,  width, color=MCOL, alpha=0.75, edgecolor='#333', linewidth=0.8, label=r'$rdh54\Delta$')
for i, pv in enumerate(pvals):
    sig_r = '***' if pv<0.001 else '**' if pv<0.01 else '*' if pv<0.05 else 'ns'
    ymax = max(w_ctrl[i], w_mut[i]) + 1.5
    ax.text(i, ymax, sig_r, ha='center', fontsize=9)
ax.set_xticks(x_pos); ax.set_xticklabels(labels_w, fontsize=9)
ax.set_xlabel('Return window after large step', fontsize=8)
ax.set_ylabel('% large steps followed\nby return (>50% back)', fontsize=8)
ax.set_title('E   Return events (Motif C-equivalent)\n"Large step → return to origin"', fontweight='bold', loc='left', fontsize=9)
ax.legend(fontsize=8); ax.tick_params(labelsize=7); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

# ── Panel F: Example trajectories annotated with large/return events ──────────
ax = fig.add_subplot(gs[2, 2:])
ax.set_facecolor('#FAFAFA')
# Pick one cell per genotype with the most return events
best_c = feat_ctrl.sort_values('return_3fr', ascending=False).iloc[0]['uid']
sub_c  = ctrl[ctrl['uid']==best_c].sort_values('frame')
xc, yc = sub_c['x_rel_px'].values, sub_c['y_rel_px'].values
dxc, dyc = np.diff(xc), np.diff(yc)
sc = np.sqrt(dxc**2 + dyc**2)

best_m = feat_mut.sort_values('return_3fr', ascending=False).iloc[0]['uid']
sub_m  = mut[mut['uid']==best_m].sort_values('frame')
xm, ym = sub_m['x_rel_px'].values, sub_m['y_rel_px'].values
dxm, dym = np.diff(xm), np.diff(ym)
sm = np.sqrt(dxm**2 + dym**2)

# Plot with offset for clarity
offset = 60
ax.plot(xc, yc, color=CCOL, lw=1.0, alpha=0.8, label='Control (best return cell)', zorder=2)
ax.plot(xm+offset, ym, color=MCOL, lw=1.0, alpha=0.8, label=r'$rdh54\Delta$ (best return cell)', zorder=2)
ax.plot(xc[0], yc[0], 'o', ms=5, color=CCOL, zorder=3)
ax.plot(xm[0]+offset, ym[0], 'o', ms=5, color=MCOL, zorder=3)

# Annotate large-step + return events (ctrl)
for i in range(len(sc)):
    if sc[i] >= THRESH_85 and i+4 < len(xc):
        dist_ret = np.sqrt((xc[i+3]-xc[i])**2 + (yc[i+3]-yc[i])**2)
        if dist_ret < 0.5 * sc[i]:  # return event
            ax.annotate('', xy=(xc[i+3], yc[i+3]), xytext=(xc[i+1], yc[i+1]),
                        arrowprops=dict(arrowstyle='<->', color='#FF8800', lw=1.5))
            ax.scatter(xc[i+1], yc[i+1], s=50, color='#FF8800', zorder=4, marker='^')

ax.set_xlabel('x (px)', fontsize=8); ax.set_ylabel('y (px)', fontsize=8)
ax.set_title('F   Example trajectories\n(arrows: large-step return events)', fontweight='bold', loc='left', fontsize=9)
ax.legend(fontsize=7.5); ax.tick_params(labelsize=7)
ax.text(0.5, -0.15, f'Threshold: {THRESH_85:.1f} px  |  pixel ≈ 0.065 μm  →  {THRESH_85*0.065:.2f} μm',
        transform=ax.transAxes, ha='center', fontsize=7.5, color='#555555')
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

# ── Panel G: RF accuracy vs permuted labels (significance test) ───────────────
ax = fig.add_subplot(gs[3, :2])
ax.set_facecolor('white')
rf_real = RandomForestClassifier(n_estimators=200, random_state=42)
real_scores = cross_val_score(rf_real, X, y, cv=skf, scoring='accuracy')
real_acc = real_scores.mean()

print("\nRunning permutation test (200 permutations)…")
perm_accs = []
for _ in range(200):
    y_perm = rng.permutation(y)
    rf_p   = RandomForestClassifier(n_estimators=100, random_state=42)
    p_sc   = cross_val_score(rf_p, X, y_perm, cv=skf, scoring='accuracy')
    perm_accs.append(p_sc.mean())
perm_accs = np.array(perm_accs)
p_perm = np.mean(perm_accs >= real_acc)
print(f"  Real RF accuracy: {real_acc:.3f} | Permuted mean: {perm_accs.mean():.3f} | p={p_perm:.4f}")

ax.hist(perm_accs, bins=25, color='#AAAAAA', alpha=0.7, edgecolor='#555555', linewidth=0.8, label='Permuted labels')
ax.axvline(real_acc, color='#CC2222', lw=2.5, label=f'Real: {real_acc:.3f} (p={p_perm:.3f})')
ax.axvline(0.5, color='#333333', lw=1, linestyle='--', label='Chance (0.5)')
ax.set_xlabel('Classification accuracy', fontsize=8)
ax.set_ylabel('Count (200 permutations)', fontsize=8)
ax.set_title('G   Permutation test (RF accuracy significance)', fontweight='bold', loc='left', fontsize=9)
ax.legend(fontsize=8); ax.tick_params(labelsize=7); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

# ── Panel H: Top features comparison (violin) ─────────────────────────────────
ax = fig.add_subplot(gs[3, 2:])
ax.set_facecolor('white')
# Top 4 discriminating features: scatter comparison
top4 = list(importances.head(4).index)
x_positions = np.arange(len(top4))
for xi, feat in enumerate(top4):
    vc = feat_ctrl[feat].values; vm = feat_mut[feat].values
    # Normalise to [0,1] for display
    vmin = min(vc.min(), vm.min()); vmax = max(vc.max(), vm.max())
    if vmax > vmin:
        vc_n = (vc - vmin) / (vmax - vmin); vm_n = (vm - vmin) / (vmax - vmin)
    else:
        vc_n = vc; vm_n = vm
    jc = rng.uniform(-0.12, 0.12, len(vc_n))
    jm = rng.uniform(-0.12, 0.12, len(vm_n))
    ax.scatter(xi - 0.2 + jc, vc_n, s=8, alpha=0.4, color=CCOL)
    ax.scatter(xi + 0.2 + jm, vm_n, s=8, alpha=0.4, color=MCOL)
    ax.plot([xi-0.2, xi-0.2], [np.percentile(vc_n,25), np.percentile(vc_n,75)], color=CCOL, lw=3, alpha=0.8)
    ax.plot([xi+0.2, xi+0.2], [np.percentile(vm_n,25), np.percentile(vm_n,75)], color=MCOL, lw=3, alpha=0.8)
    ax.plot(xi-0.2, np.median(vc_n), '_', ms=12, mew=2.5, color=CCOL)
    ax.plot(xi+0.2, np.median(vm_n), '_', ms=12, mew=2.5, color=MCOL)
    _, pv = stats.mannwhitneyu(vc, vm, alternative='two-sided')
    sig_t = '***' if pv<0.001 else '**' if pv<0.01 else '*' if pv<0.05 else 'ns'
    ax.text(xi, 1.05, sig_t, ha='center', fontsize=9)

ax.set_xticks(x_positions)
ax.set_xticklabels([t.replace('_',' ') for t in top4], fontsize=7.5, rotation=15, ha='right')
ax.set_ylabel('Normalised value', fontsize=8)
ax.set_title('H   Top 4 discriminating features\n(grey=control, red=rdh54Δ)', fontweight='bold', loc='left', fontsize=9)
ax.tick_params(labelsize=7); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

fig.suptitle('S. cerevisiae meiotic telomere trajectory ML analysis  |  Control vs rdh54Δ  |  Δt=5 min, n=46/137 cells',
             fontsize=9.5, y=0.999)

OUT = '/sessions/upbeat-awesome-wozniak/mnt/Nature com/DRAFT/cerevisiae_ML_analysis.pdf'
plt.savefig(OUT, dpi=150, bbox_inches='tight', facecolor='white')
print(f"\n✓ Figure saved → {OUT}")
plt.close()
